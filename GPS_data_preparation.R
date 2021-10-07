## GPS Data preparation 

### code to 
### 1. clean the GPS data, 
### 2. define the study area, 
### 3. remove points outside study area and when participant did not wear a GPS,
### 4. Calculate UDs and 
### 5. VI for all GPS data and night/day separated data
### 6. Calculate VI and UDs for simultaneous GPS wear data,
### 7. check week to week movement fidelity
### 8. Calculate proximity for all GPS data, coerced to one week GPS data and summarize proximity
################################################################################

source("clean_GPS_functions.R")
prj <- CRS("+proj=utm +zone=39S +ellps=WGS84 +datum=WGS84")

########################################
## 1. Clean GPS data to remove erroneous points ####
########################################
## use this on GPS data all stacked into one data frame with no points before/after participant had GPS
trimmed_stacked_dat <- read.csv(file = "trimmedGPSdata_Mandena1.csv") ## personally protected information
### social_netid = participant unique identifier
### animal_id	= NA
### gps_week	= sequential 1 = 1st week participant wore GPS, etc
### latitude	= decimal degrees as recorded by iGot-U 120
### longitude	= decimal degrees as recorded by iGot-U 120
### altitude	= meters  as recorded by iGot-U 120
### timestamp = M/d/Y H:M:S AM/PM as recorded by iGot-U 120

trimmed_stacked_dat <- trimmed_stacked_dat %>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) %>% 
  mutate(gps_week = as.character(gps_week)) %>% 
  mutate(social_netid = str_replace_all(social_netid, "-", ".")) %>% ## change A-SNH000 format to A.SNH to avoid later problems with "-"
  filter(!is.na(social_netid)) %>% 
  dplyr::select(-animal_id) %>%  #remove animal gps data 
  mutate(qi = row_number()) # give each row a unique id

step1data <- gps_flags_1(raw = trimmed_stacked_dat, 
                         timelag = 0, 
                         prj = prj, 
                         tmax = 30, 
                         smax = 30, 
                         altmin = 0, 
                         altmax = 775)
## summarize what to remove
table(step1data$flag_fri, step1data$flag_badpt)

## save it
saveRDS(step1data, "step1data.RDS")

rm(trimmed_stacked_dat)
########################################
## 2. define study area ####
########################################
library(adehabitatHR)
## using the cleaned GPS data from step 1 formated as a SpatialPointsDataFrame
step1data <- step1data %>% 
  filter(flag_badpt == FALSE) %>% #remove points flagged because of alt, loop speed, etc
  filter(flag_fri == FALSE) #remove points flagged because GPS distribution day

## make spatial
coordinates(step1data) <- ~xcoord + ycoord
proj4string(step1data) <- prj

## look at minimum convex polygon sizes
mcp.area(step1data, percent = seq(90, 100, by = 0.5)) # look for big drops in area
ct <- c(99.5, 98.5, 97, 94.5) #list of mcp percent where drops are

## look at those on a map
towns <- data.frame(town = c("Mandena", "Manantenina", "Andapa","Sarahandrano", "Marovato", "Sambava", "Marojejy Peak"),
                    latitude = c(-14.477049, -14.497213, -14.663495, -14.607567, -14.602941, -14.266516, -14.449485),
                    longitude = c( 49.814700,  49.821347,  49.651451, 49.647759, 49.634014, 50.166636,  49.732515))
towns <- st_as_sf(towns, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
towns <- st_transform(towns, prj)

plot(step1data, pch = 21, cex = 0.1) 
plot(towns, pch = 15, cex = 1, col = "red", add = T)
for (i in 1:length(ct)) {
  plot(mcp(step1data, percent = ct[i]), add = T)
}

#finer selection of polygon
mcp.area(step1data, percent = seq(94, 95, by = 0.1)) 
ct <- c(94.8, 94.7, 94.6, 94.5)
plot(mcp(step1data, percent = ct[1]))
for (i in 2:length(ct)) {
  plot(mcp(step1data, percent = ct[i]), add = T)
}
plot(step1data, pch = 21, cex = 0.1, add = T) 
plot(towns, pch = 15, cex = 1, col = "red", add = T)

## finer still
ct <- c(94.69, 94.67, 94.65, 94.63, 94.1)
plot(mcp(step1data, percent = ct[1]))
for (i in 2:length(ct)) {
  plot(mcp(step1data, percent = ct[i]), add = T)
}
plot(step1data, pch = 21, cex = 0.1, add = T) 

## Once a good contour is identified save it as the study area
studyarea <- mcp(step1data, percent = 94.63)

## plot the study area and check if a buffer is needed 
### (eg a few points in clusters on the edge of the area that are part of an "in-area" trajectory)
plot(raster::buffer(studyarea, width = 200));plot(studyarea, add = T);plot(step1data, pch = 21, cex = 0.1, add = T)

## update the study area to include the buffer if buffer is needed
studyarea <- raster::buffer(studyarea, width = 200) 
class(studyarea) #SpatialPolygons class object of the study area boundary
raster::area(studyarea)/(1000^2)

## save study area
saveRDS(studyarea, file = "gridarea_Mandena1.RDS")

rm(ct, towns)
########################################
## 3. filter points to study area and days a GPS was worn ####
########################################
step1data <- as.data.frame(step1data)
step3data <- gps_flags_2(step1data, 
                         type= "person", 
                         grid = raster::buffer(studyarea, width = 5), 
                         prj=prj) 

## summarise what was filtered
## grid:
step3data %>% 
  filter(flag_badpt == FALSE, flag_fri == FALSE) %>% 
  count(flag_grid) #584871/(584871 + 33097)
## no wear:
step3data %>% 
  filter(flag_badpt == FALSE, flag_fri == FALSE, flag_grid == FALSE) %>% 
  count(flag_day) #154675 / (154675 + 430196)

step3data <- step3data %>% 
  filter(flag_badpt == FALSE, # based on duplicates, speed, loops, altitude
         flag_fri == FALSE, # remove fridays
         flag_grid == FALSE, # remove points outside of grid
         flag_day == FALSE) %>% 
  dplyr::select(!starts_with("flag"))

## save the clean data
saveRDS(step3data, file = "CleanPersonGPS_Mandena1.RDS")

## and a list of everyone that wore a GPS
woregps <- step3data %>% 
  count(social_netid) %>%  #count number of gps fixes recorded once data cleaned
  filter(n >=240) 
write.csv(woregps, "woregps.csv", row.names = F)

## number of days of data
step3data %>% 
  filter(social_netid %in% woregps$social_netid) %>% 
  mutate(id_day = paste(social_netid, as.Date(timestamp), sep = "_")) %>% 
  pull(id_day) %>% 
  n_distinct()

## number of days of data per person
step3data %>% 
  filter(social_netid %in% woregps$social_netid) %>% 
  mutate(id_day = paste(social_netid, as.Date(timestamp), sep = "_")) %>% 
  select(social_netid, id_day) %>% 
  distinct() %>% 
  count(social_netid) %>% 
  summarise(mean_days = mean(n), median_days = median(n), sd_days = sd(n), min_days = min(n), max_days = max(n))

## number of weeks of data per person
gpswks <- data.frame(
  d1 = c("10/04/2019", "10/11/2019", "10/18/2019", "10/25/2019", "11/01/2019", "11/08/2019", "11/15/2019"),
  d2 = c("10/10/2019", "10/17/2019", "10/24/2019", "10/31/2019", "11/07/2019", "11/14/2019", "11/21/2019")
) %>% 
  mutate(d1 = as.Date(d1, format = "%m/%d/%Y", tz = "UTC"),
         d2 = as.Date(d2, format = "%m/%d/%Y", tz = "UTC"))
step3data %>% 
  mutate(dt = as.Date(timestamp)) %>% 
  select(social_netid, dt) %>% 
  distinct() %>% 
  mutate(wk = case_when(dt %in% seq(gpswks$d1[1], gpswks$d2[1], by = 1) ~ 1,
                        dt %in% seq(gpswks$d1[2], gpswks$d2[2], by = 1) ~ 2,
                        dt %in% seq(gpswks$d1[3], gpswks$d2[3], by = 1) ~ 3,
                        dt %in% seq(gpswks$d1[4], gpswks$d2[4], by = 1) ~ 4,
                        dt %in% seq(gpswks$d1[5], gpswks$d2[5], by = 1) ~ 5,
                        dt %in% seq(gpswks$d1[6], gpswks$d2[6], by = 1) ~ 6,
                        dt %in% seq(gpswks$d1[7], gpswks$d2[7], by = 1) ~ 7)) %>% 
  select(-dt) %>% 
  distinct() %>% 
  count(social_netid) %>% 
  summarise(mean_wks = mean(n), median_wks = median(n), sd_wks = sd(n), min_wks = min(n), max_wks = max(n))

########################################
## 4. Calculate entire UDs  ####
########################################

## packages
library(tidyverse)
library(move)
library(maptools)
library(sp)

## parameters for dbbmm
prj <- CRS("+proj=utm +zone=39S +ellps=WGS84 +datum=WGS84") #set the CRS projection
loc.er = 9.2 ## the accuracy of the GPS in meters
## margin and window size are the number of fixes (eg mar = 5 means 13 minutes)
mar = 5 # The margin used for the behavioral change point analysis. 
win.siz = 19 # The size of the moving window along the track.This number has to be odd.
ext = 0.5 # extension of the bbox of the track, using default
ts = 3/15 # time interval of each fix in min, best to set so random freq fixes don't increase comp time
bbgrid <- NA #the bounding box argument (doing this so if want can use code below for the grid)

########################################
### 4a. prepare the data  ####
########################################
## load data from step3 
step3data <- readRDS(file = "CleanPersonGPS_Mandena1.RDS") %>% dplyr::select(-geometry)

step3data <- step3data  %>%   
  dplyr::select(social_netid, timestamp, xcoord, ycoord)
## make sure data in order (cleanest form of GPS point data from step 3)
step3data <- step3data[order(step3data$social_netid, step3data$timestamp),]

## grid building using the studyarea
grid.poly <- buffer(studyarea, width = 5000) #add a huge buffer because gaps result in huge 100% UDs
x <- seq(round(grid.poly@bbox[1]), round(grid.poly@bbox[3]), by=10) #warning that 1m object is huge! probably a much better way to do this
y <- seq(round(grid.poly@bbox[2]), round(grid.poly@bbox[4]), by=10)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
raster::crs(xy) <- prj
gridded(xy) <- TRUE
class(xy) ## SpatialPixels

## save those for VI calc
saveRDS(xy, "studyareagrid_Mandena1_10m.RDS")

## make into a raster for dbbmm and ud calc
cells <- raster::raster(xy, layer=0, values=FALSE)

rm(studyarea, grid.poly, x, y)

########################################
### 4b. dBBMM of UD ####
########################################

## to check if working use a small subset of step3data called "df"
# df <- step3data[step3data$social_netid %in% c("ID1", "ID2"),] ##chose a couple of individuals to replace ID1, ID2

## once happy with parameters run on the whole thing
## rename step3data as df 
df <- step3data

## make into a move object
df.mov <- move(x=df$xcoord, y=df$ycoord, time=df$timestamp, proj=prj, animal=df$social_netid)
## can do more data cleaning to a move object if you want (see vignette)

## to deal with the gappiness of the tracks (which results in a ud that is a blob)
## calculate the variance of the gappy track
var <- brownian.motion.variance.dyn(df.mov, location.error = loc.er, 
                                    margin = mar, window.size = win.siz) #takes a long time! good to save the object so don't have to rerun

## tell it to ignore all segments with a timelag greater than 15 min AND all segments that are islands (segment before and after ignored)
var@interest[unlist(lapply(timeLag(var, units = "mins"), c, NA))>30] <-FALSE #tell it to 'ignore' segments where >15min between fixes
var@interest[lead(var@interest == FALSE) & lag(var@interest  == FALSE)] <- FALSE # remove standalone points (those with a removed point before and after)

## calculate the dbbmm
memory.size(memory.limit()*3) #increase memory alloted to R so don't get Error: cannot allocate vector
dbbmm <- brownian.bridge.dyn(var, location.error = loc.er, margin = mar, window.size = win.siz,
                             ext = ext, time.step = ts, raster = cells, bbox = bbgrid)
## calculate the UD
ud <- getVolumeUD(dbbmm)

########################################
### 4c. plot 95%UD ####
########################################
## so just see the 95% UD colored
ud95 <- ud
ud95[ud95>0.95]<- NA
## to see 1 individual's plot
namelist <- names(dbbmm)
who <- "A.SNH006" #chose an individual on the namelist
lay <- which(namelist == who)
plot(ud95[[lay]], main = who) #ud95
# plot(ud[[lay]], main = who) #whole ud
plot(var[var@trackId == who]@data[,c("x", "y")], col = "red", cex = 0.05, pch = 20)
contour(ud[[lay]], levels = c(.5, .95, 0.99), add = T, lwd = 0.05, lty = c(2, 1, 3))

########################################
### 4d. ## save both dbbmm and ud  ####
########################################

## save as RDS:
saveRDS(dbbmm, file = "dbbmm_Mandena1.RDS")
saveRDS(ud, file = "ud_Mandena1.RDS")

## save ESRI shapefile (if want to see in arc)
# library(rgdal)
# library(raster)
# for (i in 1:length(names(dbbmm))) {
#   ## the rasters (ArcMap can make the contours on its own)
#   writeRaster(ud[[i]], filename = paste("myfilepath", names(dbbmm)[i], ".tif", sep = ""), overwrite=TRUE)
#   # con1 <- raster2contour(dbbmm[[i]], level = c(.25, .3, .35, .4, .45, .50, .55, .6, .65, .7, .75, .8, .85, .90, .95))
#   # writeOGR(con1, dsn = "", layer = names(dbbmm[i]), driver = "ESRI Shapefile", overwrite_layer = TRUE)
# }

########################################
## 5. VI of UDs ####
########################################
## need these objects
xy <- readRDS("studyareagrid_Mandena1_10m.RDS")
dbbmm <- readRDS("dbbmm_Mandena1.RDS")
ud <- readRDS("ud_Mandena1.RDS")

########################################
### 5a. VI function ####
########################################

## this is adapted from adehabitatHR to work with a dBBMM
## homerange/utiliz-distribution/vi_calc.R

my.kerneloverlap.m1 <- function(dbbmm, dbbmm.ud, grid, 
                                percent=95, conditional=FALSE, ...)
  #dbbmm is the output from the brownian.bridge.dyn lapply mod1 (list of rasters)
  #dbbmm.ud is the output from method 1 mod1.ud
  #grid is the study area selected, must be gridded not raster (xy)
  #see adehabitatHR:kerneloverlapHR for other parameters and overlap calculations
{
  ## Verifications
  # method <- match.arg(method) #select method to use
  if (ncol(coordinates(xy))>2) #check just xy in coords
    stop("xy should be defined in two dimensions")
  if (percent >100) #check percent
    stop("percent should not be >100")
  if (percent <=1) {
    stop("percent should be >1") 
  }
  
  ## set up
  x <-dbbmm
  vol <- dbbmm.ud
  gp <- gridparameters(grid) #make sure not integer(0)
  
  ids <- names(x)
  res <- data.frame()
  for (i in 1:(length(ids)-1)) {
    for (j in (i+1):length(ids)) {
      tmp <- data.frame(id1 = ids[i], id2 = ids[j])
      res <-rbind(res, tmp)
    }
  }
  
  ## loop through each pair
  for (q in 1:nrow(res)) { #if use directional method need to change these to 1:length(x) for i and j
    i = res$id1[q]
    j = res$id2[q]
    
    vi <- values(x[[i]])
    vj <- values(x[[j]])
    ai <- values(vol[[i]])*100 
    aj <- values(vol[[j]])*100
    ai[ai<=percent] <- 1
    ai[ai>percent] <- 0
    aj[aj<=percent] <- 1
    aj[aj>percent] <- 0
    
    if (conditional) {
      vi <- vi*ai
      vj <- vj*aj
      r = sum(pmin(vi, vj))*(gp[1,2]^2)
      res[q,3] <- r
    } 
    else {
      r = sum(pmin(vi, vj))*(gp[1,2]^2)
      res[q,3] <- r
    }
    
  } 
  
  colnames(res) <- c("id1", "id2", paste("VI", percent, sep = "_"))
  return(res)
} 

########################################
### 5b. VI of entire UDs ####
########################################

vi <- my.kerneloverlap.m1(dbbmm = dbbmm, dbbmm.ud = ud, grid = xy, percent = 95, conditional = TRUE)
vi50 <- my.kerneloverlap.m1(dbbmm = dbbmm, dbbmm.ud = ud, grid = xy, percent = 50, conditional = TRUE)

viboth <- vi %>% full_join(vi50)
summary(viboth)

## save the file for use in models later on
write.csv(viboth, file = "vi_all_Mandena1.csv", row.names = F) 

rm(dbbmm, ud, vi, vi50, viboth)
########################################
### 5c. VI of separate night and day UDs ####
########################################
## if need to reload data, also need my.kerneloverlap.m1 function from(5a)
step3data <- readRDS(file = "CleanPersonGPS_Mandena1.RDS")  %>% dplyr::select(-geometry)
step3data <- step3data %>% 
  dplyr::select(social_netid, timestamp, xcoord, ycoord)
## make sure data in order (cleanest form of GPS point data from step 3)
step3data <- step3data[order(step3data$social_netid, step3data$timestamp),]

## repeat the grid cells from above
xy <- readRDS("studyareagrid_Mandena1_10m.RDS")
## make into a raster for dbbmm and ud calc
cells <- raster::raster(xy, layer=0, values=FALSE)

## parameters for dbbmm (same as above)
prj <- CRS("+proj=utm +zone=39S +ellps=WGS84 +datum=WGS84") #set the CRS projection
loc.er = 9.2 ## the accuracy of the GPS in meters
## margin and window size are the number of fixes (eg mar = 5 means 13 minutes)
mar = 5 # The margin used for the behavioral change point analysis. 
win.siz = 19 # The size of the moving window along the track.This number has to be odd.
ext = 0.5 # extension of the bbox of the track, using default
ts = 3/15 # time interval of each fix in min, best to set so random freq fixes don't increase comp time
bbgrid <- NA #the bounding box argument (doing this so if want can use code below for the grid)

## create weekday and time of day columns to subset the data with
## note that you can also do this for days of the week and hours of the day just change names in loop
main <- step3data %>% 
  mutate(timestamp.tmp = timestamp) %>% 
  tidyr::separate(timestamp.tmp, into = c("orgDate", "orgHr"), sep = " ") %>% 
  mutate(wDay = lubridate::wday(orgDate, label = TRUE)) %>% 
  mutate(wDay = as.character(wDay)) %>% 
  tidyr::separate(orgHr, into = c("hr", "min", NA), sep = ":") %>% 
  mutate(daynight = ifelse(hr %in% c("20", "21", "22", "23", "00", "01", "02", "03", "04", "05"), "night", "day")) %>% ##using Oct dark hours roughly...
  dplyr::select(social_netid, timestamp, xcoord, ycoord, hr, min, wDay, daynight)

## subset the data
tod <- unique(main$daynight)

## create new df for storage
vi.tod <- data.frame(id1= character(), id2 = character())

for (i in 1:length(tod)) { 
  toi <- tod[i]
  df1 <- main %>% filter(daynight == toi)
  df.mov <- move(x=df1$xcoord, y=df1$ycoord, time=df1$timestamp, proj=prj, animal=df1$social_netid)
  df.mov <- df.mov[which(df.mov@trackId %in% names(which(n.locs(df.mov) > win.siz)))]# remove traj with fewer fixes than the window.size
  var <- brownian.motion.variance.dyn(df.mov, location.error = loc.er, margin = mar, window.size = win.siz)
  
  var@interest[unlist(lapply(timeLag(var, units = "mins"), c, NA))>15] <-FALSE #tell it to 'ignore' segments where >15min between fixes
  var@interest[lead(var@interest == FALSE) & lag(var@interest  == FALSE)] <- FALSE # remove standalone points (those with a removed point before and after)
  
  dbbmm <- brownian.bridge.dyn(var, location.error = loc.er, margin = mar, window.size = win.siz,
                               ext = ext, time.step = ts, raster = cells, bbox = bbox(cells))
  ud <- getVolumeUD(dbbmm)
  # saveRDS(dbbmm, paste0(toi, "_Mandena1.RDS"))
  
  vi95 <-my.kerneloverlap.m1(dbbmm, ud, xy, percent=95, conditional = TRUE)
  vi50 <-my.kerneloverlap.m1(dbbmm, ud, xy, percent=50, conditional = TRUE)
  vi <- full_join(vi95, vi50)
  colnames(vi) <- c("id1", "id2", paste("VI_95", toi, sep="_"), paste("VI_50", toi, sep = "_"))
  
  vi.tod <- vi.tod %>% full_join(vi, by = c("id1", "id2"))
  
  
  rm(df1, df.mov, var, dbbmm, ud, vi, vi95, vi50)
}
summary(vi.tod) 

## save the file for use in models later on
write.csv(vi.tod, file = "vi_nightday_Mandena1.csv", row.names = F)

########################################
## 6. VI of GPS cowear UDs ####
########################################
## for using GLM need VI of time both people were wear a GPS simultaneously
## using the same parameters as above to calculate UD and VI (see step 4)
## homerange/utiliz-distribution/dBBMM-overlap_cowearGPS.R

## make a list of all possible id combos
woregps <- woregps[,1]
pairs <- data.frame()
for (i in 1:(length(woregps)-1)) {
  for (j in (i+1):length(woregps)) {
    tmp <- data.frame(id1 = woregps[i], id2 = woregps[j])
    pairs <-rbind(pairs, tmp)
  }
}


########################################
### 6a. write cowear VI function ####
########################################

cowear <- function(x, pairids, locs, xy, min.overlap, prj, loc.er, mar, win.siz, ext, ts, cells, bbgrid){ 
  
  ## make a separate dataframe for each individual in each pair
  print(x)
  vi50 = c()
  vi95 = c()
  id1 <- pairids[,1]
  id2 <- pairids[,2]
  
  df1 <- locs[locs$social_netid == id1[x],]
  df1$t <- difftime(df1$timestamp, as.POSIXct("2019-01-01 00:00:00", tz = "UTC"), units = "secs") ## make time in to seconds from Jan 1, 2019 00:00:00 for fuzzyjoin
  
  df2 <- locs[locs$social_netid == id2[x],]
  df2$t <- difftime(df2$timestamp, as.POSIXct("2019-01-01 00:00:00", tz = "UTC"), units = "secs") ## make time in to seconds from Jan 1, 2019 00:00:00 for fuzzyjoin
  
  ## join the dataframes based on "fuzzy" matched times
  close_times <- 
    fuzzyjoin::difference_inner_join(df1, df2, by = "t", max_dist = 60*60) #max distance of 1 hour
  
  ## make sure have a minimual number of observations for each person
  p1 = close_times %>% dplyr::select(ends_with("x")) %>% distinct() %>% nrow()
  p2 = close_times %>% dplyr::select(ends_with("y")) %>% distinct() %>% nrow()
  
  if (p1 > min.overlap & p2 > min.overlap) {
    ## separate back to own df's
    df3 <- close_times %>% 
      dplyr::select(ends_with("x"), -t.x) %>% 
      rename(social_netid = social_netid.x, timestamp = timestamp.x, xcoord = xcoord.x, ycoord = ycoord.x)
    df3 <- close_times %>% 
      dplyr::select(ends_with("y"), -t.y)%>% 
      rename(social_netid = social_netid.y, timestamp = timestamp.y, xcoord = xcoord.y, ycoord = ycoord.y) %>% 
      rbind(df3) %>% 
      distinct()
    
    ## calc UD for each indiv
    df.mov <- move(x=df3$xcoord, y=df3$ycoord, time=df3$timestamp, proj=prj, animal=df3$social_netid)
    ## calculate the variance of the gappy track
    var <- brownian.motion.variance.dyn(df.mov, location.error = loc.er, margin = mar, window.size = win.siz) 
    var@interest[unlist(lapply(timeLag(var, units = "mins"), c, NA))>15] <-FALSE #tell it to 'ignore' segments where >15min between fixes
    var@interest[lead(var@interest == FALSE) & lag(var@interest  == FALSE)] <- FALSE # remove standalone points (those with a removed point before and after)
    
    segs <- as.data.frame(var) %>% dplyr::count(trackId, interest) %>% dplyr::filter(interest)
    if (segs[1,3] > min.overlap & segs[2,3] > min.overlap) { #need at least 12 windows
      dbbmm <- brownian.bridge.dyn(var, location.error = loc.er, margin = mar, window.size = win.siz,
                                   ext = ext, time.step = ts, raster = cells, bbox = bbgrid)
      ## calculate the UD
      ud <- getVolumeUD(dbbmm)
      
      ## calc VI
      v1 <-my.kerneloverlap.m1(dbbmm, ud, xy, method="VI", percent=50, conditional = TRUE)[1,3]
      
      v2 <-my.kerneloverlap.m1(dbbmm, ud, xy, method="VI", percent=95, conditional = TRUE)[1,3]
      
      ## save variables for df
      vi50 = c(vi50, v1)
      vi95 = c(vi95, v2)
      
      rm(df3, df.mov, ls.mov, var, mod1, mod1.ud, vi50, vi95)
      
    }else{
      v1 = NA
      v2 = NA
      vi50 = c(vi50, v1)
      vi95 = c(vi95, v2)
    }
  } else {
    v1 = NA
    v2 = NA
    vi50 = c(vi50, v1)
    vi95 = c(vi95, v2)
  }
  
  results = data.frame(id1 = id1[x], id2 = id2[x], VI_50 = v1, VI_95 = v2)
  return(results)
}

########################################
### 6b. Calc cowear VI ####
########################################
## using the same data and parameters reloaded at the start of step 5c.
## this parameter was different
ext = 0.3 # extension of the bbox of the track was different, 
#### ext is not considered when the raster is a (non-numeric) raster in brownian.bridge.dyn()

## WARNING THIS TOOK 45.25 HOURS TO RUN 
## I recommend running in subsets or might try to make a cluster, mclapply, but didn't seem to make faster with my trials
memory.size(memory.limit()*3) # to stop error cannot allocate vector of size X
system.time(
  {vi <- lapply(1:nrow(pairs), cowear, pairids = pairs, locs = step3data, xy = xy, min.overlap = 240, 
                prj, loc.er, mar, win.siz, ext, ts, cells, bbgrid)
  vi1 <- as.data.frame(do.call(rbind, vi))})

## look at summary 
summary(vi1)

## some checks
vi1 %>% n_distinct() #good all are distinct!
length(which(is.na(vi1$VI_95)))/nrow(vi1) ## half of the pairs are NAs

## save it!
write.csv(vi1, file = "vi_cowear_Mandena1.csv", row.names = F)

########################################
## 7. Movement fidelity week to week ####
########################################
## load data again
step3data <- readRDS(file = "CleanPersonGPS_Mandena1.RDS")  %>% dplyr::select(-geometry)
## make sure data in order (cleanest form of GPS point data from step 3)
step3data <- step3data[order(step3data$social_netid, step3data$timestamp),]

## also using those same parameters except cells because less costly if just set to 10m size instead of use raster
ext = 0.3
cells = 10

########################################
### 7a. figure out who wore a GPS of >3 weeks and >10 days ####
########################################

df.summary <- step3data %>% 
  group_by(social_netid) %>% 
  summarise(nweeks = n_distinct(gps_week), ndays = n_distinct(lubridate::yday(timestamp))) %>% 
  arrange(desc(nweeks)) %>% 
  filter(nweeks>3, ndays >10)
df.summary
nrow(df.summary)

df.subset <-  step3data %>% 
  filter(social_netid %in% df.summary$social_netid)

rm(step3data, df.summary)

########################################
### 7b. calculate weekly UDs for those people  ####
########################################

##dbbmm
df.mov <- move(x=df.subset$xcoord, y=df.subset$ycoord, time=df.subset$timestamp, proj=prj, animal=df.subset$id_week)

ls.mov <- split(df.mov)
var <- brownian.motion.variance.dyn(df.mov, location.error = loc.er, 
                                    margin = mar, window.size = win.siz) #takes a long time! good to save the object so don't have to rerun

## tell it to ignore all segments with a timelag greater than 15 min AND all segments that are islands (segment before and after ignored)
var@interest[unlist(lapply(timeLag(var, units = "mins"), c, NA))>15] <-FALSE #tell it to 'ignore' segments where >15min between fixes
var@interest[lead(var@interest == FALSE) & lag(var@interest  == FALSE)] <- FALSE # remove standalone points (those with a removed point before and after)

## calculate the dbbmm
memory.size(memory.limit()*3) #increase memory alloted to R so don't get Error: cannot allocate vector
dbbmm <- brownian.bridge.dyn(var, location.error = loc.er, margin = mar, window.size = win.siz,
                             ext = ext, time.step = ts, raster = cells, bbox = bbgrid)
## calculate the UD
ud <- getVolumeUD(dbbmm)

########################################
### 7c. week to week VI  ####
########################################
## VI_95
vi <-my.kerneloverlap.m1(dbbmm = dbbmm, dbbmm.ud = ud, grid = xy, percent = 95, conditional = TRUE)

vi.self <- vi %>% 
  filter(id1 != id2) %>% #remove error line from function (no overlap with self in same week)
  separate(col = id1, into = c("id1", "sweek"), sep = "_") %>% 
  separate(col = id2, into = c("id2", "rweek"), sep = "_") %>% 
  filter(id1 == id2) %>% #get only overlaps with self (different weeks)
  dplyr::select(-id2)

vi.summary <- vi.self %>% 
  group_by(id1) %>% 
  summarise(avg.vi = mean(VI_95), nweeks = n())
vi.summary
mean(vi.summary$avg.vi);sd(vi.summary$avg.vi);range(vi.summary$avg.vi)

## VI_50
vi.core <-my.kerneloverlap.m1(dbbmm = dbbmm, dbbmm.ud = ud, grid = xy, percent = 50, conditional = TRUE)

vicore.self <- vi.core %>% 
  filter(sender != receiver) %>% #remove error line from function (no overlap with self)
  separate(col = sender, into = c("sender", "sweek"), sep = "_") %>% 
  separate(col = receiver, into = c("receiver", "rweek"), sep = "_") %>% 
  filter(sender == receiver) %>% #get only overlaps with self (different weeks)
  dplyr::select(-receiver)

vicore.summary <- vicore.self %>% 
  group_by(sender) %>% 
  summarise(avg.vi = mean(VI_50), nweeks = n())
vicore.summary
mean(vicore.summary$avg.vi);sd(vicore.summary(avg.vi));range(vicore.summary$avg.vi)

########################################
## 8. Proximity ####
########################################
## Close-contacts/proximity.R ####
########################################
### 8a. Find distance threshold ####
########################################
## calculate Spatial threshold value (dc)
SpThValues<-contact::findDistThresh(n = 100000, acc.Dist1 = loc.er, ## LE of GPS = 9.2+/- .2m
                                    pWithin1 = 95, spTh = 0) #spTh represents the initially-defined spatial threshold for contact. Note that we've chosen to use 100,000 in-contact point-location pairs here.
SpThValues  #Note that because these confidence intervals are obtained from distributions generated from random samples, every time this function is run, results will be slightly different.
dc <- unname(SpThValues[["spTh.adjustments"]][2])
dc #pTh_98%Capture number =17.03
rm(SpThValues)

########################################
### 8b. proximity function ####
########################################
## function of calculate simultaneous, proximal (within distance threshold set above),
### proximity/prox = n proximal / n sim
### also can calculated time spent together and if 'point' or 'continual' contacts
##### this last part is commented out because it takes ~40 hours to run
myprox <- function(pts, dc, fix_freq, window_min){
  
  require(fuzzyjoin)
  if (window_min < fix_freq) {
    stop("window_min must be >= to fix_freq")
  }
  
  tc = fix_freq/2 #the temporal threshold must be 1/2 of the fix frequency or less for this to work
  
  ## be careful that dataframe has the columns in this order regardless of column name
  colnames(pts) <- c("id", "xcoord", "ycoord", "timestamp")
  
  ## make time in to seconds from Jan 1, 2019 00:00:00
  pts <- pts %>% 
    dplyr::mutate(timestamp = difftime(timestamp, as.POSIXct("2019-01-01 00:00:00", tz = "UTC"), units = "secs"))
  
  idslist <- unique(pts$id)
  ids <- data.frame()
  for (i in 1:(length(idslist)-1)) {
    for (j in (i+1):length(idslist)) {
      tmp <- data.frame(id1 = idslist[i], id2 = idslist[j])
      ids <-rbind(ids, tmp)
    }
  }
  
  out <- data.frame()
  
  for (i in 1:nrow(ids)) {
    ## make a separate dataframe for each individual in each pair
    df1 <- pts %>% 
      dplyr::filter(id == ids$id1[i])
    df2 <- pts %>% 
      dplyr::filter(id == ids$id2[i])
    
    ## join the dataframes based on "fuzzy" matched times (timestamp within +/- tc)
    close_times <- 
      fuzzyjoin::difference_full_join(df1, df2, by = "timestamp", max_dist = tc) 
    
    ## calculations to get prox and other info out
    if(nrow(na.omit(close_times)) > 1){ ## only do the calculations if there is temporal overlap
      close_times <- close_times %>% 
        dplyr::arrange(timestamp.x) %>% #full_join puts NAs at the bottom, want in temporal order
        dplyr::mutate(d= sqrt(((xcoord.x-xcoord.y)^2)+((ycoord.x-ycoord.y)^2))) %>% #calc distance between points
        dplyr::mutate(spa = dplyr::if_else(d <=dc, 1, 0)) #%>% #check if distance in <= distance threshold (dc)
      # dplyr::mutate(con = dplyr::if_else(spa & lag(spa) & abs(lag(timestamp.x) - timestamp.x) <= (fix_freq + 60), 1,0)) #is interaction a single point or many? #ADDED IN SOME ERROR FOR THE FIX_FREQ HERE
      
      ###### Add in to calculate elapsed time ######
      # this basically creates a vector to count 'triggers' (consecutive camera clicks or GPS fixes)
      # and another vector to count 'sequences' that are subtriggers within each trigger
      
      # trigger2 <- data.frame() # empty dfs
      # sequence_n <- data.frame()
      # x <- 0
      # seq_n <- 1
      
      # cloogy loop
      # for(j in 1:length(close_times$con)) {
      #   if ( is.na(close_times$con[j]) ) { # dealing with NAs
      #     x = x+1
      #     trig=x
      #     seq_n = 1
      #   }
      #   else {
      #     if (close_times$con[j] == 0) {
      #       x = x+1
      #       trig = x
      #       seq_n = 1
      #     }
      #     else {
      #       trig = x
      #       seq_n = seq_n+1
      #     }
      #   }
      #   trigger2 = rbind(trigger2, trig)
      #   sequence_n = rbind(sequence_n, seq_n)
      # }
      
      
      # close_times$Trigger = trigger2$X1 # add the output to the df
      # close_times$Seq = sequence_n$X1
      
      ## calculate the elapsed time for each trigger
      # elapsed = close_times %>% 
      #   filter(spa ==1) %>% 
      #   dplyr::mutate(mean_time = (as.numeric(timestamp.x) + as.numeric(timestamp.y))/2)%>% #avg start and end time of the trigger
      #   dplyr::group_by(Trigger) %>% 
      #   dplyr::summarise(elapsed = last(mean_time)-first(mean_time)) %>% 
      #   dplyr::mutate(elapsed = elapsed + tc) ## assume all contacts last for a minimum of tc (so point contacts don't last for 0 min)
      # elapsed2 = as.numeric(sum(elapsed$elapsed))/60 #total time they are together in min
      # if (nrow(elapsed)>0) {
      #   avgelapsed = elapsed2/length(na.omit(unique(elapsed$Trigger)))
      # }else{ avgelapsed=0}
      
      # elapsed_thresh = elapsed %>% 
      # dplyr::filter(elapsed >= window_min) ## to get only elapsed times above threshold
      # elapsed_thresh =  as.numeric(sum(elapsed_thresh$elapsed))/60 #total time above threshold (might want to change this to number of unique triggers above threshold)
      
      nsim = nrow(na.omit(close_times)) #length(na.exclude(close_times$mean_time))
      nspa = sum(close_times$spa, na.rm = TRUE) #number of spatially proximal and temporally simultaneous fixes based on thresholds
      # ncon = sum(close_times$con, na.rm = TRUE) #number of fixes where the previous fix was also spatially proximal and temporally simultaneous
      prox = nspa/nsim #proportion of spatially proximal and temporally simultaneous fixes over all temporally simultaneous fixes
      # mprox = nspa*fix_freq/60 #the number of minutes of overlap (assuming fix frequency was constant)
      
      tmp <- data.frame(id1 = unique(df1$id), id2 = unique(df2$id),
                        n_simultaneous = nsim, n_proximal = nspa, 
                        # n_continual_contact = ncon, avg_contact_time = avgelapsed, #proximal_minutes = mprox,
                        # sim_pt_first = sim_1st, sim_pt_last = sim_last, 
                        # prox_mins_above_threshold = elapsed_thresh,#elap,
                        prox = prox)
      out <- rbind(out, tmp)
    }
    
    else{ #return NAs for pairs without temporal overlap
      nsim = 0
      nspa = NA
      # ncon = NA
      prox = NA
      
      # avgelapsed = NA
      
      tmp <- data.frame(id1 = unique(df1$id), id2 = unique(df2$id),
                        n_simultaneous = nsim, n_proximal = nspa, 
                        # n_continual_contact = ncon, 
                        # avg_contact_time = avgelapsed, 
                        # prox_mins_above_threshold = NA, 
                        prox = prox)
      out <- rbind(out, tmp)
    }
    
  }
  return(out)
}

########################################
## 8c. Find observed proximity ####
########################################
memory.size(memory.limit()*3)
step3data <- readRDS(file = "CleanPersonGPS_Mandena1.RDS")
step3data <- step3data %>% 
  dplyr::select(social_netid, xcoord, ycoord, timestamp) ## arrange the dataframe to only have needed columns in correct order

## set the fix frequency (3 min) in seconds, this can be upped and still works eg any one who over laps within 10 min instead of 3
## the temporal threshold is half the fix freq - 3 min is 1.5 min on either side of the fix
fix_freq = 3*60

## set the window min (this is the threshold in 'total time above threshold' in output) to the fix frequency as well
window_min = 3*60

## run the function (takes  ~1hr)
prox.df <- myprox(pts = step3data, dc = dc, fix_freq = fix_freq, window_min = window_min)

head(prox.df)
summary(prox.df)

## save it!
write.csv(prox.df, file = "prox_3min_Mandena1.csv", row.names = F)

########################################
### 8d. Find 'pooled' proximity ####
########################################
## change the dates so that all trajectories occur over the same week
## this assumes people have a routine behavior week to week (see 7)
## this will only capture interactions that occur on the same weekday (eg all Mondays)
## so there will be some NAs if people don't have matching days of wear (eg person 1 has only Saturday and person 2 has only Sunday)
## but people with that minimal amount of overlap don't seem very valuable to model builiding anyways
## one way to get rid of these people is filter out by overlap (say need >1/2 day or >240 fixes)

cdat <- step3data %>% 
  tidyr::separate(timestamp, into = c("orgDate", "orgHr"), sep = " ") %>% #separate date and time
  mutate(wDay = lubridate::wday(orgDate, label = TRUE)) %>% #get the weekday of each date
  mutate(wDay = as.character(wDay)) %>% 
  mutate(newDate = case_when(wDay == "Sat" ~ "2019-01-05", #make all the same weekdays the same date
                             wDay == "Sun" ~ "2019-01-06",
                             wDay == "Mon" ~ "2019-01-07",
                             wDay == "Tue" ~ "2019-01-08",
                             wDay == "Wed" ~ "2019-01-09",
                             wDay == "Thu" ~ "2019-01-10",)) %>% 
  unite(timestamp, c(newDate, orgHr), sep = " ") %>% #join the newdate to the origional time
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) %>% #make into time format again
  dplyr::select(social_netid, xcoord, ycoord, timestamp) #select the relevant columns

## using the same parameters as above
# dc = 17.03
# fix_freq = 3*60
# window_min = 3*60

## re run the prox function
### this took ~2.5 hrs
pooled.prox <- myprox(pts = cdat, dc = dc, fix_freq = fix_freq, window_min = window_min)

## breifly compare pooled and observed prox
summary(pooled.prox)
summary(prox.df)

## save it!
write.csv(pooled.prox, file = "coercedprox_3min_Mandena1.csv", row.names = F)

########################################
### 8e. Summarise proximity output ####
########################################
## Individuals that wore a GPS at the same time 
prox.1.summary <- prox.df %>% 
  pivot_longer(cols = c(-id1, -id2), names_to = "column", values_to = "val") %>% 
  group_by(column) %>% 
  summarise(min = min(val, na.rm=TRUE), max = max(val, na.rm=TRUE), median = median(val, na.rm=TRUE), mean = mean(val, na.rm=TRUE), sd = sd(val, na.rm=TRUE))
prox.1.summary

## Individuals that cowore a GPS for at least half a day (>=240 pts)   
prox.240.summary <- prox.df %>% 
  filter(n_simultaneous >=240)%>% 
  pivot_longer(cols = c(-id1, -id2), names_to = "column", values_to = "val") %>% 
  group_by(column) %>% 
  summarise(min = min(val, na.rm=TRUE), max = max(val, na.rm=TRUE), median = median(val, na.rm=TRUE), mean = mean(val, na.rm=TRUE), sd = sd(val, na.rm=TRUE))
prox.240.summary
