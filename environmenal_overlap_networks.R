## Environmental overlap networks  

########################################
## 1. Calculate use of areas  ####
########################################
require(tidyverse)
require(rgdal)
require(raster)
require(maptools)
require(sf)
require(exactextractr)

woregps <- read.csv("woregps.csv")[,1]

## people's utilization distributions
dbbmm <- readRDS(file = "dbbmm_Mandena1.RDS") #from GPS_data_preparation 4d.
dbbmm <- dbbmm[[which(names(dbbmm) %in% woregps)]]
 
## classified raster (done in ArcMap)
r <- raster("../Corrected_10m_SVM/Corrected_v2/Reclassified_2021_SVM_FUL.tif")

## check extent and CRS of both
dbbmm@crs
r@crs
## reproject ud95 to same crs as r (reproject ud95 because continuous, if reproject r need to use method = "ngb" because categorical and will induce NAs)
prj <- r@crs
dbbmm <- projectRaster(dbbmm, crs = prj)
dbbmm@crs
r@crs

########################################
## 1.1 95% UD ####
########################################
## subset to the 95% use area (the 100% is nonsense probability area from irregularities in GPS fix rate)
ud95 <- raster()
for(i in 1:nlayers(dbbmm)){
  d <- dbbmm[[i]] #extract one layer
  d[d == 0] <- NA #make the cells with prop UD of 0 into NA
  vec = na.omit(values(d)) #create vector of all non-NA values
  b = cumsum(vec[order(vec, decreasing = T)]) #order that vector and sum along it to find value cutoff for 95% sum
  line = length(b[b<=0.95])
  # line = length(b[b<=0.05]) #get the ordered location of where that cutoff is
  thresh = vec[order(vec, decreasing = T)][line] #find the value based on the location
  d[d < thresh] <- NA #set all cells less than that value to NA
  # sum(values(d), na.rm = T)
  ud95 <- stack(ud95, d) #save each layer to a raster stack
}

nlayers(ud95)
rm(d, b, i, line, thresh, vec, dbbmm)

########################################
## 1.2 Calculate use and cover type in each cell ####
########################################

## clip ud95 and r to study area using polygon from GPS_data_preparation step 2.
grid.poly <- readRDS(file = "gridarea_Mandena1.RDS")
grid.poly <- spTransform(grid.poly, prj) #convert to the same crs as ud95 and r

## plot to show stuff
image(r)
plot(log10(ud95[[1]]), add = T)
plot(grid.poly, add = T)

# Redo with the grid shifted 10m in x and y directions  
# This should cover any affects due to the arbitrary grid placement

## first get the extent of the original grid.poly
egp <- round(as.vector(extent(grid.poly))) #add + or - # here to shift grid slightly!

## make a list of all possible shifts:
ax <- c(0, 10, 10, -10, -10, 10, -10, 0, 0)
ay <- c(0, 10, -10, 10, -10, 0, 0, 10, -10)

pplcoverlist <- list()

for (i in 1:length(ax)) {
  dx = ax[i]
  dy = ay[i]
  
  xmin = egp[1] + dx
  xmax = egp[2] + dx
  ymin = egp[3] + dy
  ymax = egp[4] + dy
  ## get the difference between 30 and the remainder when divide by 30
  addx = 30 - (xmax - xmin)%%30
  addy = 30- (ymax - ymin)%%30
  ## add the difference to the max directions
  xmax = xmax + addx
  ymax = ymax + addy
  ## make a new (square) polygon to crop r and ud95 to
  x = c(xmin, xmax, xmax, xmin, xmin)
  y = c(ymin, ymin, ymax, ymax, ymin)
  new.grid.poly <- Polygon(cbind(x, y))
  new.grid.poly <- Polygons(list(new.grid.poly), 1)
  new.grid.poly <- SpatialPolygons(list(new.grid.poly), proj4string = prj)
  
  ##removed the check divisible by 30 step
  rm(xmin, xmax, ymin, ymax, addx, addy, x, y)
  
  #clip the raster and ud to that polygon (can also clip the ud95)
  rclip <- crop(r, new.grid.poly)
  # plot(rclip); plot(new.grid.poly, lwd = 2, add = T)
  
  #make the sampling grid
  grid <- st_make_grid(rclip, cellsize = 30) 
  grid <- st_sf(grid)
  
  ##calculate UD in each cell
  use <-  exact_extract(ud95, grid, 'sum', append_cols=TRUE) # and use ud clip here
  names(use) = gsub(pattern = "sum.", replacement = "", x=names(use))
  
  ## check column sums to about 0.95
  use %>% 
    dplyr::select(starts_with("A")) %>% 
    summarise(across(.cols = everything(), .fns = sum)) %>% 
    pivot_longer(cols = everything()) %>% 
    summarise(range= range(value)) %>% 
    print()
  
  ## calc land class pixels in that cell
  coversum <- list()
  for (j in 1:7) {
    r0 <- rclip+1
    r0[r0 != j] <- NA
    r0[r0 == j] <- 1
    coversum[j+1] <- exactextractr::exact_extract(r0, grid, 'sum', append_cols=TRUE)
  }
  coversum <- as.data.frame(do.call(cbind, coversum))
  colnames(coversum) <- c("Primary", "Secondary", "Rice", "Village", "Other", "Water", "Bare_ground")
  
  ## join use and coversum
  pplcover <- cbind(use, coversum)
  ## make the grid cell id a column
  pplcover <- pplcover %>% rownames_to_column(var = "gridcell")
  
  ## save that to a list
  pplcoverlist[[i]] <- pplcover
}

## save all of those
saveRDS(pplcoverlist, file = "landclass_30m_overlap_Mandena1.RDS")

########################################
## 1.3 Useful checks ####
########################################
## Look at a small section of the raster and grid overlayed to make sure 
### each grid cell has 9 raster cells. Ignore the edges of the raster that is 
### due to the extent argument I used to clip the raster  
# plot(grid);plot(rclip, add = T)
cent = round(centroid(new.grid.poly))
ext = c(cent[1]-90, cent[1]+90, cent[2]-90, cent[2]+90)
plot(rclip, ext = ext);plot(grid, extent = 
                              st_bbox(c(xmin =ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4]), crs = prj), 
                            add = T)
rm(ext, cent, new.grid.poly)

# basic summaries
## Look at the summaries of everything and max UD proportion values by landcover type and use  

## the total areas covered by each land class in 10m-sq
pplcoverlist[[1]] %>% 
  pivot_longer(cols = c(Primary, Secondary, Rice, Village, Other, Water, Bare_ground), names_to = "covertype", values_to = "coverprop") %>% 
  group_by(covertype) %>% 
  summarise(classified_area_10msq = sum(coverprop))

pplcoverlist[[1]] %>% 
  mutate(gridcell_use = rowSums(across(starts_with("A"), na.rm = T))) %>% 
  dplyr::select(gridcell, gridcell_use, Primary, Secondary, Rice, Village, Other, Water, Bare_ground) %>% 
  arrange(desc(gridcell_use)) %>% head(n = 10L)

## check column sums to about 0.95 (also in the loop as a check)
pplcoverlist[[1]] %>% 
  dplyr::select(starts_with("A")) %>% 
  summarise(across(.cols = everything(), .fns = sum)) %>% 
  pivot_longer(cols = everything()) %>% 
  summarise(range= range(value))

########################################
## 1.4 Plot use by cover type ####
########################################
w = pplcoverlist[[1]] %>% dplyr::select(gridcell, Water)

w$xcoor = rep(seq(1:311), 310) #x dimension of grid esp[2]-esp[1]
w$ycoor = rep(1:310, each=311) #y dimension of grid esp[4]-esp[3]

ggplot(w, aes(x=xcoor, y=ycoor))+
  geom_tile(aes(fill=Water))

rm(w)
########################################
## 2. Make networks ####
########################################
library(tidyverse)
library(igraph)

pplcoverlist <- readRDS(file = "landclass_30m_overlap_Mandena1.RDS") #from end of step 1.2
## take the average value for each grid cell
### this takes forever - probably because it is close to 1GB
# olap = plyr::aaply(plyr::laply(pplcoverlist, as.matrix), c(2, 3), mean)

## this is MUCH faster
pplcovermatrices <- list()
gridcell <- pplcoverlist[[1]]$gridcell
for (i in 1:length(pplcoverlist)) {
  pplcovermatrices[[i]] <- as.matrix(pplcoverlist[[i]][,-1])
  rownames(pplcovermatrices[[i]]) <- gridcell
}

olap <- Reduce('+', pplcovermatrices)/length(pplcovermatrices)

olap <-  as.data.frame(olap)
olap <-  olap %>% rownames_to_column("gridcell")

rm(pplcoverlist, pplcovermatrices, gridcell)

## check this by plotting (like above)
w = olap %>% select(gridcell, Water)
w$xcoor = rep(seq(1:311), 310) #x dimension of grid esp[2]-esp[1]
w$ycoor = rep(1:310, each=311) #y dimension of grid esp[4]-esp[3]
ggplot(w, aes(x=xcoor, y=ycoor))+
  geom_tile(aes(fill=Water))
rm(w)

########################################
## 2.1 write functions ####
########################################

## function to remove grid cells that have no human use and make table into long format
shape_raster = function(df){
  
  # calculate rowsums
  df$total = df %>% select(2:ncol(df)) %>% rowSums()
  
  df=df %>%
    filter(total>0) %>% # retain nonzeros
    select(gridcell, contains("SNH")) %>% 
    pivot_longer(cols=contains("SNH")) %>% #SNH is in the uid for all people
    return()
}

# two functions -- one for village, one for nonvillage
find_dominant_cover = function(wide_df, long_df){
  # match to gridcell composition 
  olap_l_f_j = long_df %>%
    left_join(select(wide_df, gridcell, Primary, Secondary, Rice, Village, Other, Water, Bare_ground),
                         by="gridcell")
  
  # Determine dominant cover type
  domcov = olap_l_f_j %>%
    pivot_longer(cols=c(Primary,Secondary, Rice, Village, Other, Water, Bare_ground), 
                 names_to="Cover", values_to="CellCt") %>%
    group_by(gridcell, name) %>% 
    filter(CellCt==max(CellCt)) %>%
    ungroup()
  
  # filter out rows that don't show any use ()
  domcov_f = domcov %>% filter(value > 0)
  w = data.frame(gridcell = wide_df$gridcell,
                 xcoor = rep(seq(1:310), 311), #cast back to grid dimensions
                 ycoor = rep(1:311, each=310))
  domcov_f = domcov_f %>% 
    left_join(w, by = "gridcell")
  return(domcov_f)
}
# nonvillage 
find_dominant_cover_nv = function(wide_df, long_df){
  # match to gridcell composition
  olap_l_f_j = long_df %>%
    left_join(select(wide_df, gridcell, Primary, Secondary, Rice, Village, Other, Water, Bare_ground),
                         by="gridcell") %>%
    filter(Village==0)# remove any cells with village

  # Determine dominant cover type

  domcov = olap_l_f_j %>%
    pivot_longer(cols=c(Primary,Secondary, Rice, Village, Other, Water, Bare_ground), names_to="Cover", values_to="CellCt") %>%
    group_by(gridcell, name) %>%
    filter(CellCt==max(CellCt)) %>%
    ungroup()

  domcov_f = domcov %>% filter(value > 0)
  w = data.frame(gridcell = wide_df$gridcell,
                 xcoor = rep(seq(1:310), 311),
                 ycoor = rep(1:311, each=310))

  domvoc_f = domcov_f %>%
    left_join(w, by = "gridcell")
  return(domcov_f)
}


up_land = function(gph, threshpct = 0.9){
  g_bp_w = bipartite.projection(gph)
  g_human_w = g_bp_w$proj2
  V(g_human_w)%>%length();E(g_human_w)%>%length() 
  
  E(g_human_w)$width = log(E(g_human_w)$weight)
  thresh = quantile(E(g_human_w)$weight, threshpct)
  
  E(g_human_w)$wt1 = E(g_human_w)$width %>%as.data.frame()%>% mutate(m=ifelse(.>log(threshpct),0,.)) %>%select(m)
  E(g_human_w)$wt2 = E(g_human_w)$width %>%as.data.frame()%>% mutate(m=ifelse(.<log(threshpct),0,.)) %>%select(m)
  
  return(g_human_w)
}

## for specific land cover types weight edges by prop cover and use in gridcell
landtype_viz_weighted = function(df, covertype, colcol){
  
  # create column weighting UD by cell count proportion
  df$value2 = df$value * df$CellCt/9
  
  # now filter for any positive UD for the given landcover type
  w_olap = df %>% filter(Cover == covertype & value2 > 0) 
  
  g_w = w_olap %>% 
    select(gridcell, name, value2) %>% 
    graph_from_data_frame(., directed = F)
  
  V(g_w)$type = bipartite_mapping(g_w)$type 
  E(g_w)$weight = E(g_w)$value2
  
  summary(V(g_w)$type)
  land_cells = as.numeric(c(summary(V(g_w)$type)[2]))
  people = as.numeric(c(summary(V(g_w)$type)[3]))
  length(E(g_w)$weight)
  
  V(g_w)$label = c(rep("", land_cells), V(g_w)$name[(land_cells+1):(land_cells+people)])  
  
  
  return(g_w)
}

########################################
## 2.2 form networks ####
########################################
# Several individuals lived in areas surrounded by rice. The landcover classification missed
# a couple of these. We added a nominal amount of 'village' to allow for potential filtering:
olap %>%
  select(A.SNH039, gridcell,Rice, Village, Secondary, Primary, Bare_ground, Other) %>% 
  filter(A.SNH039 >0.05) %>% 
  arrange(desc(A.SNH039)) 

# if more than 5% of UD in a single 30m2 cell, give a small village weight
## visual inspection showed these individuals had a small house inside of rice cells so we reclassified their house as village manually
olap[which(olap$A.SNH039 > 0.05 & olap$Village==0),]$Village = 0.1
olap[which(olap$A.SNH088 > 0.05 & olap$Village==0),]$Village = 0.1
olap[which(olap$A.SNH072 > 0.05 & olap$Village==0),]$Village = 0.1

# One individual seemed to have a high UD around where we *think* is camp
# exclude camp gridcells, but need to double-check these with coordinates!
# gridcells_excluded = olap[which(olap$A.SNH072 > 0.02 & olap$Rice>0),]$gridcell

## apply the functions
olap_long = shape_raster(olap) 
domcov_f = find_dominant_cover(wide_df = olap, long_df = olap_long) 
domcov_f_nv = find_dominant_cover_nv(wide_df = olap, long_df = olap_long)

## make a list of all the humanIDs
humanIDs = olap %>% select(contains("SNH")) %>% names()

######## make bipartite graph ######################
# all_bp_nz = all_bp(long_df = olap_long)
all_bp_nz <- olap_long %>% 
  filter(value > 0) %>% 
  select(gridcell, name, value) %>% 
  graph_from_data_frame(., directed = F)
## tell igraph it is bipartite
V(all_bp_nz)$type = bipartite_mapping(all_bp_nz)$type 
## values are edge weights
E(all_bp_nz)$weight = E(all_bp_nz)$value
land_cells = as.numeric(c(summary(V(all_bp_nz)$type)[2]))
people = as.numeric(c(summary(V(all_bp_nz)$type)[3]))
V(all_bp_nz)$label = c(rep("", land_cells), V(all_bp_nz)$name[(land_cells+1):(land_cells+people)])  
rm(land_cells, people)

###### and its unipartite projection ################
all_uni_nz = up_land(gph = all_bp_nz)

###### repeat for each land cover type #######
## for the MS we only looked at rice
#### rice
rice_bp = landtype_viz_weighted(df = domcov_f_nv, covertype = "Rice", colcol = "dodgerblue")
rice_uni = up_land(rice_bp, threshpct = .9)

########################################
## 2.3 Weighted unipartite projections ####
########################################
####### All land cover types ############
# all bp -- note that this is nonzero
all_proj = bipartite_projection(all_bp_nz)

# retrieve edgelist
V(all_proj$proj2)$name
comparisons = combn(V(all_proj$proj2)$name,2, simplify=T) %>% t()

edgedf = as_data_frame(all_bp_nz, "edges")
df = edgedf %>% 
  select(from, to, value) %>% 
  pivot_wider(names_from = "from", values_from = "value", values_fill=0) # takes a second

df = dplyr::rename(df, name2=to)
df_num = df[,-1]

head(comparisons)
dim(comparisons)

# create empty vector
vec=c()

## product of proportion cover and proportion ud for each gridcell (about 10 mins)
extent = length(comparisons)/2
for(i in 1:extent){
  v = apply(df_num[df$name2%in%(comparisons[i,]),], MARGIN = 2, FUN=prod) %>% sum()
  vec=c(vec,v)
}

all_uni_nz
checklist=paste(comparisons[,1],comparisons[,2], sep="|")
which( (checklist == attr(E(all_uni_nz),"vnames")) == F) # should be 0
checklist[match(attr(E(all_uni_nz),"vnames"),checklist)] == attr(E(all_uni_nz),"vnames")
E(all_uni_nz)$weight = vec[match(attr(E(all_uni_nz),"vnames"),checklist)]

plot(all_uni_nz, edge.width = E(all_uni_nz)$weight*100, layout=layout_with_mds, vertex.label.cex = 0.7)

############# Rice ###########
rice_proj = bipartite_projection(rice_bp)

# retrieve edgelist
V(rice_proj$proj2)$name
comparisons = combn(V(rice_proj$proj2)$name,2, simplify=T) %>% t()

edgedf = as_data_frame(rice_bp, "edges")
df = edgedf %>% 
  select(from, to, value2) %>% 
  pivot_wider(names_from = "from", values_from = "value2", values_fill=0) 

df = dplyr::rename(df, name2=to)
df_num = df[,-1]

head(comparisons)
dim(comparisons)

# create empty vector
vec=c()

## product of proportion cover and proportion ud for each gridcell (~2 mins)
extent = length(comparisons)/2
for(i in 1:extent){
  v = apply(df_num[df$name2%in%(comparisons[i,]),], MARGIN = 2, FUN=prod) %>% sum()
  vec=c(vec,v)
}

# need to 
rice_uni
checklist=paste(comparisons[,1],comparisons[,2], sep="|")
which( (checklist == attr(E(rice_uni),"vnames")) == F) # should be 0 if arranged the same way
checklist[match(attr(E(rice_uni),"vnames"),checklist)] == attr(E(rice_uni),"vnames") # all true
E(rice_uni)$weight = vec[match(attr(E(rice_uni),"vnames"),checklist)]

plot(rice_uni, edge.width = E(rice_uni)$weight*100, layout=layout_with_mds, vertex.label.cex = 0.7)

####### Save the networks ###########
saveRDS(list(enviro_all = all_uni_nz, 
             enviro_rice = rice_uni), file = "networks_enviro.RDS") 

########################################
## 3. network wide summaries ####
########################################
## write function to calculate
widestats = function(graph, name, type = c("uni", "bi")){
  ## need matrix in this format for bipartite
  m <- as.matrix(igraph::as_adjacency_matrix(graph, type = "both", attr = "weight"))
  m <- m[which(rownames(m) %in% humanIDs), -which(colnames(m) %in% humanIDs)] 
  widestatsdf = data.frame(
    network=  name,
    connected_strong = igraph::is.connected(graph, mode = "strong"),
    connected_weak = igraph::is.connected(graph, mode = "weak"),
    diameter = igraph::diameter(graph),
    avg_distance = igraph::average.path.length(graph),
    density = igraph::graph.density(graph),
    transitivity_global = igraph::transitivity(graph, type = "global"))
  if(type == "uni"){
    louvain_modularity = igraph::modularity(igraph::cluster_louvain(graph))
    widestatsdf$louvain_modularity = louvain_modularity
  }
  if(type == "bi"){
    bipartite_modularity = bipartite::computeModules(m, forceLPA = TRUE)@likelihood #takes awhile
    widestatsdf$biparite_modularity = bipartite_modularity
  }
  
  return(widestatsdf)
}

## bipartite network summary 
ALL = widestats(all_bp_nz, "ALL", "bi")
RICE = widestats(rice_bp, "RICE", "bi")

(graphmetrics = rbind(ALL, RICE)) 

write.csv(graphmetrics, "widestats_enviro_bp.csv", row.names=F)

## unipartite network summary
## fix issue with missing nodes (missing due to subsetting method):
vcount(all_uni_nz); vcount(rice_uni)

#rice
nn <-  setdiff(V(all_uni_nz)$name, V(rice_uni)$name)
rice1 <- add.vertices(rice_uni, nv = length(nn), name = nn)

ALL = widestats(all_uni_nz, "ALL", "uni")
RICE = widestats(rice_uni, "RICE", "uni")

graphmetrics = rbind(ALL, RICE)

## add in number of disconnected nodes:
graphmetrics$n_disconnected <- 123 - c(vcount(all_uni_nz), vcount(rice_uni))
graphmetrics

## save it
write.csv(graphmetrics, "widestats_enviro_uni.csv", row.names=F)

########################################
## 4. Unipartite centrality scores ####
########################################
############# All land cover types #############
e_cent_a = igraph::eigen_centrality(all_uni_nz)$vector
b_cent_a = igraph::betweenness(all_uni_nz, weights = (1 - E(all_uni_nz)$weight), directed=F)
d_cent_a = igraph::degree(all_uni_nz)
s_cent_a = igraph::strength(all_uni_nz)

e_cent_a = data.frame(e_cent_a) %>%
  mutate(humanID = rownames(.))
b_cent_a = data.frame(b_cent_a) %>%
  mutate(humanID = rownames(.))
d_cent_a = data.frame(d_cent_a) %>% 
  mutate(humanID = rownames(.))
s_cent_a = data.frame(s_cent_a) %>% 
  mutate(humanID = rownames(.))

cent_a = e_cent_a %>% 
  left_join(b_cent_a) %>% 
  left_join(d_cent_a) %>% 
  left_join(s_cent_a)

names(cent_a) = c("EIGEN","humanID","BETWEEN","DEGREE","STRENGTH")
cent_a$Network = "ALL"

cent_a_long = cent_a %>% 
  pivot_longer(cols=c(EIGEN,BETWEEN:STRENGTH))
cent_a_long

############# Rice #############
e_cent_r = igraph::eigen_centrality(rice_uni)$vector
b_cent_r = igraph::betweenness(rice_uni, weights = (1 - E(rice_uni)$weight), directed=F)
d_cent_r = igraph::degree(rice_uni)
s_cent_r = igraph::strength(rice_uni)

e_cent_r = data.frame(e_cent_r) %>%
  mutate(humanID = rownames(.))
b_cent_r = data.frame(b_cent_r) %>%
  mutate(humanID = rownames(.))
d_cent_r = data.frame(d_cent_r) %>% 
  mutate(humanID = rownames(.))
s_cent_r = data.frame(s_cent_r) %>% 
  mutate(humanID = rownames(.))

cent_r = e_cent_r %>% 
  left_join(b_cent_r) %>% 
  left_join(d_cent_r) %>% 
  left_join(s_cent_r)

names(cent_r) = c("EIGEN","humanID","BETWEEN","DEGREE","STRENGTH")
cent_r$Network = "RICE"

cent_r_long = cent_r %>% 
  pivot_longer(cols=c(EIGEN,BETWEEN:STRENGTH))
cent_r_long

########################################
## 4.1 combine unipartite centrality scores ####
########################################
all_centrality = rbind(cent_a_long, cent_r_long)
head(all_centrality)
names(all_centrality) = c("ID","network","metric","value")

all_centrality %>%
  group_by(metric, network) %>% 
  mutate(value= scale(value)) %>% 
  ungroup %>% 
  ggplot(aes(x=reorder(ID, -value), y=network))+
  geom_tile(aes(fill=kader:::cuberoot(value)))+
  facet_wrap(~metric, nrow=5)+
  scale_fill_gradient2(low=scales::muted("blue"), mid="white", high=scales::muted("red"))+
  theme(aspect.ratio= 1/20, axis.text.x = element_text(angle=90, hjust=1, vjust=0))

all_centrality %>%
  group_by(metric, network) %>% 
  mutate(value= scale(value)) %>% 
  ungroup %>% 
  ggplot(aes(x=reorder(ID, -value), y=metric))+
  geom_tile(aes(fill=kader:::cuberoot(value)))+
  facet_wrap(~network, nrow=4)+
  scale_fill_gradient2(low=scales::muted("blue"), mid="white", high=scales::muted("red"))+
  theme(aspect.ratio= 1/20, axis.text.x = element_text(angle=90, hjust=1, vjust=0))

## save it
write.csv(all_centrality, "centrality_enviro.csv", row.names = F)

