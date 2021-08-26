## Compare networks

## plot, edge weight vs naming, correlation in centrality
################################################################################

library(tidyverse)
library(igraph)
woregps <- read.csv("woregps.csv")[,1]
########################################
## 1. Plot networks ####
########################################

## to save all plots use
# pdf(file = "", bg = "transparent") 
# #### OR
# png(file = "", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
# plot
# dev.off()

set.seed(34727)
## Import all the networks  ####
nd_net <- readRDS("network_name_nondirected.RDS")
d_net <- readRDS("network_name_directed.RDS")
enviro_net <- readRDS("networks_enviro.RDS")

cc.mean <- read.csv("close_contact_edges_nonames.csv") %>% 
  select(id1, id2, imputed.weight)
cc_net_nonames <- graph_from_data_frame(cc.mean, vertices = woregps)

## set the color palettes ####
col.notnamed <- adjustcolor("#1B9E77", alpha.f = 0.5) #greenish
col.named <- "#D95F02" #orange
col.recip <- "#7570B3" #purple
col.vert <-  adjustcolor("#666666", alpha.f = 0.85) #grey
col.edges <- adjustcolor("#1B9E77", alpha.f = 0.5) #also greenish

## full directed naming network ####
### layout to use for all graphs
l = layout_with_fr(d_net$n.named)

## change colors to match boxplots and other networks
V(d_net$n.named)$color <- col.vert
E(d_net$n.named)[which(E(d_net$n.named)$mutual == "grey50")]$color <-  col.named 
E(d_net$n.named)[which(E(d_net$n.named)$mutual == "blue")]$color <-  col.recip

## colored by naming relation:
# pdf(file = "plots/naming_directed_recipcolored.pdf", bg = "transparent") 
# png(file = "plots/naming_directed_recipcolored.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(d_net$n.named,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.arrow.size = 0.2,
     edge.width = E(d_net$n.named)$weight,
     edge.color = E(d_net$n.named)$color)#adjustcolor("gray", alpha.f = 0.5))
# dev.off()

## colored green edges
# pdf(file = "plots/naming_directed_green.pdf", bg = "transparent") 
# png(file = "plots/naming_directed_green.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(d_net$n.named,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.arrow.size = 0.2,
     edge.width = E(d_net$n.named)$weight,
     edge.color = col.edges)
# dev.off()

## save as df to use on other networks
nn.col <- as_data_frame(d_net$n.named, "edges")[,c("from", "to", "color")]


## close contact network ####
### for plotting purposes (otherwise network is too dense to see anything)
#### remove edges with less than 20 contacts per perfect week, round to 100th place

cc_net_plot <- cc.mean %>% 
  mutate(weight = ifelse(imputed.weight <= 20/(20*24*7), 0, round(imputed.weight, 2))) %>%
  filter(weight > 0) %>%  # remove edges that = 0
  left_join(nn.col, by = c("id1" = "from", "id2" = "to")) %>% #join to full directed naming network to get colors
  mutate(color = replace_na(color, col.notnamed)) %>% #the NAs in color did not name eachother so give them a color
  arrange(color) #so named edges are on top and easier to see

## make back into a graph
cc_net_plot <- graph_from_data_frame(cc_net_plot, vertices = woregps, directed = FALSE)
is.simple(cc_net_plot) 

## change verticies color match naming network
V(cc_net_plot)$color <- col.vert

## colored by naming edges
# pdf(file = "plots/closecontact_directed_nonames_recipcolored.pdf", bg = "transparent")
# png(file = "plots/closecontact_directed_nonames_recipcolored.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(cc_net_plot,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = E(cc_net_plot)$weight*5,
     edge.color = E(cc_net_plot)$color)
# dev.off()

## colored green edges
# pdf(file = "plots/closecontact_directed_nonames_green.pdf", bg = "transparent")
# png(file = "plots/closecontact_directed_nonames_green.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(cc_net_plot,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = E(cc_net_plot)$new_cc*5,
     edge.color = col.edges) 
# dev.off()

rm(cc_net_plot)

## enviro all network ####

envall <- as_data_frame(enviro_net$enviro_all, "edges")
envall <- envall %>% 
  select(from, to, weight) %>% 
  mutate(weight = ifelse(weight <= 0.01^2, 0, weight)) %>% #at least 1% of both UDs overlapping
  filter(weight > 0) %>% 
  left_join(nn.col, by = c("from", "to"))  %>% #get the colors
  mutate(color = replace_na(color, col.notnamed)) %>% 
  arrange(color) # sort edges so edges in naming network are on top

#make back into a network
envall <- graph_from_data_frame(envall, vertices = woregps, directed = F) 

## change verticies color match naming network
V(envall)$color <- col.vert

## colored by naming edges
# pdf(file = "plots/enviroAll_directed_recipcolored.pdf", bg = "transparent")
# png(file = "plots/enviroAll_directed_recipcolored.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(envall,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = log10(1+E(envall)$weight)*100,
     edge.color = adjustcolor(E(envall)$color, alpha.f = 0.75))
# dev.off()

## colored green edges
# pdf(file = "plots/enviroAll_directed_green.pdf", bg = "transparent")
# png(file = "plots/enviroAll_directed_green.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(envall,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = log10(1+E(envall)$weight)*100,
     edge.color = adjustcolor(col.edges, alpha.f = 0.25))
# dev.off()

rm(envall)
## enviro rice network ####
envrice <- as_data_frame(enviro_net$enviro_rice, "edges")
envrice <- envrice %>% 
  select(from, to, weight) %>% 
  mutate(weight = ifelse(weight <= 0.001^2, 0, weight)) %>% #at least 0.1% of both UDs overlapping
  filter(weight > 0) %>% 
  left_join(nn.col, by = c("from", "to"))  %>% #get the colors
  mutate(color = replace_na(color, col.notnamed)) %>% 
  arrange(color) # sort edges so edges in naming network are on top

#make back into a network
envrice <- graph_from_data_frame(envrice, vertices = woregps, directed = F) 

## change verticies color match naming network
V(envrice)$color <- col.vert

##colored by naming edges
# pdf(file = "plots/enviroRice_directed_recipcolored.pdf", bg = "transparent")
# png(file = "plots/enviroRice_directed_recipcolored.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(envrice,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = log10(1+E(envrice)$weight)*100,
     edge.color = E(envrice)$color)
# dev.off()

## colored green edges
# pdf(file = "plots/enviroRice_directed_green.pdf", bg = "transparent")
# png(file = "plots/enviroRice_directed_green.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
plot(envrice,
     layout = l,
     vertex.size = 7,
     vertex.label = NA,
     vertex.frame.color = "#969696",#slightly lighter grey
     edge.width = log10(1+E(envrice)$weight)*100,
     edge.color = col.edges)
# dev.off()
rm(envrice)

########################################
## 2. label reciprocal edges ####
########################################
## Check mutual, which_mutual, and dyad_cencus match
ecount(d_net$n.named)
sum(which_mutual(d_net$n.named))/2 #recip
sum(!which_mutual(d_net$n.named))
dyad_census(d_net$n.named)
as_data_frame(d_net$n.named, "edges") %>% 
  mutate(named = ifelse(mutual == "blue", 2, 1)) %>% 
  select(from, to, named) %>% 
  count(named)
132/2

## use the big edge table to get n.named and recip.named
edgedf <- read.csv("big_edgelist_Mandena1.csv") %>% 
  select(id1, id2, n.named, recip.name) %>% 
  filter(id1 %in% woregps, id2 %in% woregps)
table(edgedf$recip.name)

## recip.named from freetime network 
E(d_net$freetime)$recip.name = which_mutual(d_net$freetime)
recip <- as_data_frame(d_net$freetime) %>% 
  select(from, to, recip.name) %>% 
  mutate(freetime = ifelse(recip.name, 1, 0)) %>%
  select(from, to, freetime)

edgedf <- edgedf %>% 
  left_join(recip, by = c("id1" = "from", "id2" = "to")) %>% 
  mutate(freetime = ifelse(n.named > 0 & is.na(freetime), 0, freetime),
         freetime = case_when(freetime == 0 ~ "named",
                                         freetime == 1 ~ "recip.named",
                                         TRUE ~ "not.named"))

## recip.named from farm network
E(d_net$farm)$recip.name = which_mutual(d_net$farm)
recip <- as_data_frame(d_net$farm) %>% 
  select(from, to, recip.name) %>% 
  mutate(farm = ifelse(recip.name, 1, 0)) %>%
  select(from, to, farm)

edgedf <- edgedf %>% 
  left_join(recip, by = c("id1" = "from", "id2" = "to")) %>% 
  mutate(farm = ifelse(n.named > 0 & is.na(farm), 0, farm),
         farm = case_when(farm == 0 ~ "named",
                                         farm == 1 ~ "recip.named",
                                         TRUE ~ "not.named"))

## recip.named from food network
E(d_net$food)$recip.name = which_mutual(d_net$food)
recip <- as_data_frame(d_net$food) %>% 
  select(from, to, recip.name) %>% 
  mutate(food = ifelse(recip.name, 1, 0)) %>%
  select(from, to, food)

edgedf <- edgedf %>% 
  left_join(recip, by = c("id1" = "from", "id2" = "to")) %>% 
  mutate(food = ifelse(n.named > 0 & is.na(food), 0, food),
         food = case_when(food == 0 ~ "named",
                                         food == 1 ~ "recip.named",
                                         TRUE ~ "not.named"))

########################################
## 3. df of edge weights ####
########################################

edgedf <- cc.mean %>% 
  filter(imputed.weight > 0) %>% 
  rename(cc_weight = imputed.weight) %>% 
  full_join(edgedf, by = c("id1", "id2")) 

## enivro all edge weights (unipartite)
edgedf <- as_data_frame(enviro_net$enviro_all, "edges") %>% 
  select(from, to, weight) %>% 
  rename("allE" = "weight", id1 = from, id2 = to) %>% 
  right_join(edgedf, by = c("id1", "id2")) 

## enivro rice edge weights (unipartite)
edgedf <- as_data_frame(enviro_net$enviro_rice, "edges") %>% 
  select(from, to, weight) %>% 
  rename("riceE" = "weight", id1 = from, id2 = to) %>% 
  right_join(edgedf, by = c("id1", "id2")) 

range(edgedf$cc_weight, na.rm = T)
range(edgedf$allE, na.rm = T)
range(edgedf$riceE, na.rm = T)

########################################
## 4. Edge weights and recip.named ####
########################################
##Look at number of dyads in each subgroup
table(edgedf$recip.name)
table(edgedf$freetime)
table(edgedf$farm)
table(edgedf$food)

## a list of the pairwise comparisons to make in box plots below
mycomparisons <- list(c("not.named", "named"), c("named", "recip.named"), c("not.named", "recip.named"))

## short function to make boxplots
## box plots
my_box_funt <- function(net, naming, my.title){
  print(
    ggplot(edgedf, aes(x = {{naming}}, y = {{net}}, fill = {{naming}})) + 
      geom_boxplot() + 
      # ggpubr::stat_compare_means(method = "kruskal.test") + #might need to add label.x and label.y if use this 
      ggpubr::stat_compare_means(method = "wilcox.test", comparisons = mycomparisons) +
      # ggpubr::stat_compare_means(label = "p.signif", method = "t.test", ref.group = "recip.named") + 
      theme_minimal() +
      labs(title = my.title) +
      scale_fill_brewer(palette = "Dark2")
  )
}


########################################
## 4a. close contact comparisons ####
########################################
edgedf <- edgedf %>% 
  mutate(across(c(recip.name, freetime, farm, food), ~factor(.x, levels = c("not.named", "named", "recip.named"))))
## close contact vs all named 
my_box_funt(net = cc_weight, 
            naming = recip.name, 
            my.title  = "close-contact no names - complete naming")
wilcox.test(edgedf$cc_weight[which(edgedf$recip.name == "not.named")], edgedf$cc_weight[which(edgedf$recip.name == "named")])
wilcox.test(edgedf$cc_weight[which(edgedf$recip.name == "not.named")], edgedf$cc_weight[which(edgedf$recip.name == "recip.named")])
wilcox.test(edgedf$cc_weight[which(edgedf$recip.name == "named")], edgedf$cc_weight[which(edgedf$recip.name == "recip.named")])

## close contact vs freetime 
my_box_funt(net = cc_weight, 
            naming = freetime, 
            my.title  = "close-contact no names - freetime naming")
wilcox.test(edgedf$cc_weight[which(edgedf$freetime == "not.named")], edgedf$cc_weight[which(edgedf$freetime == "named")])
wilcox.test(edgedf$cc_weight[which(edgedf$freetime == "not.named")], edgedf$cc_weight[which(edgedf$freetime == "recip.named")])
wilcox.test(edgedf$cc_weight[which(edgedf$freetime == "named")], edgedf$cc_weight[which(edgedf$freetime == "recip.named")])

## close contact vs farm 
my_box_funt(net = cc_weight, 
            naming = farm, 
            my.title  = "close-contact no names - farm naming")
wilcox.test(edgedf$cc_weight[which(edgedf$farm == "not.named")], edgedf$cc_weight[which(edgedf$farm == "named")])
wilcox.test(edgedf$cc_weight[which(edgedf$farm == "not.named")], edgedf$cc_weight[which(edgedf$farm == "recip.named")])
wilcox.test(edgedf$cc_weight[which(edgedf$farm == "named")], edgedf$cc_weight[which(edgedf$farm == "recip.named")])

## close contact vs food 
my_box_funt(net = cc_weight, 
            naming = food, 
            my.title  = "close-contact names - food naming")
wilcox.test(edgedf$cc_weight[which(edgedf$food == "not.named")], edgedf$cc_weight[which(edgedf$food == "named")])
wilcox.test(edgedf$cc_weight[which(edgedf$food == "not.named")], edgedf$cc_weight[which(edgedf$food == "recip.named")])
wilcox.test(edgedf$cc_weight[which(edgedf$food == "named")], edgedf$cc_weight[which(edgedf$food == "recip.named")])

## table of mean and sd weights
edgedf %>% 
  group_by(recip.name) %>% 
  summarise(mean_wt = mean(cc_weight, na.rm = T), sd_wt = sd(cc_weight, na.rm = T))
edgedf %>% 
  group_by(freetime) %>% 
  summarise(mean_wt = mean(cc_weight, na.rm = T), sd_wt = sd(cc_weight, na.rm = T))
edgedf %>% 
  group_by(farm) %>% 
  summarise(mean_wt = mean(cc_weight, na.rm = T), sd_wt = sd(cc_weight, na.rm = T))
edgedf %>% 
  group_by(food) %>% 
  summarise(mean_wt = mean(cc_weight, na.rm = T), sd_wt = sd(cc_weight, na.rm = T))

########################################
## 4b. enviro all comparisons ####
########################################
## enviro all vs all named 
my_box_funt(net = allE, 
            naming = recip.name, 
            my.title  = "enviro all - complete naming")
wilcox.test(edgedf$allE[which(edgedf$recip.name == "not.named")], edgedf$allE[which(edgedf$recip.name == "named")])
wilcox.test(edgedf$allE[which(edgedf$recip.name == "not.named")], edgedf$allE[which(edgedf$recip.name == "recip.named")])
wilcox.test(edgedf$allE[which(edgedf$recip.name == "named")], edgedf$allE[which(edgedf$recip.name == "recip.named")])

## enviro all vs freetime 
my_box_funt(net = allE, 
            naming = freetime, 
            my.title  = "enviro all - freetime naming")
wilcox.test(edgedf$allE[which(edgedf$freetime == "not.named")], edgedf$allE[which(edgedf$freetime == "named")])
wilcox.test(edgedf$allE[which(edgedf$freetime == "not.named")], edgedf$allE[which(edgedf$freetime == "recip.named")])
wilcox.test(edgedf$allE[which(edgedf$freetime == "named")], edgedf$allE[which(edgedf$freetime == "recip.named")])

## enviro all vs farm 
my_box_funt(net = allE, 
            naming = farm, 
            my.title  = "enviro all - farm naming")
wilcox.test(edgedf$allE[which(edgedf$farm == "not.named")], edgedf$allE[which(edgedf$farm == "named")])
wilcox.test(edgedf$allE[which(edgedf$farm == "not.named")], edgedf$allE[which(edgedf$farm == "recip.named")])
wilcox.test(edgedf$allE[which(edgedf$farm == "named")], edgedf$allE[which(edgedf$farm == "recip.named")])

## enviro all vs food 
my_box_funt(net = allE, 
            naming = food, 
            my.title  = "enviro all - food naming")
wilcox.test(edgedf$allE[which(edgedf$food == "not.named")], edgedf$allE[which(edgedf$food == "named")])
wilcox.test(edgedf$allE[which(edgedf$food == "not.named")], edgedf$allE[which(edgedf$food == "recip.named")])
wilcox.test(edgedf$allE[which(edgedf$food == "named")], edgedf$allE[which(edgedf$food == "recip.named")])

## table of mean and sd weights
edgedf %>% 
  group_by(recip.name) %>% 
  summarise(mean_wt = mean(allE, na.rm = T), sd_wt = sd(allE, na.rm = T))
edgedf %>% 
  group_by(freetime) %>% 
  summarise(mean_wt = mean(allE, na.rm = T), sd_wt = sd(allE, na.rm = T))
edgedf %>% 
  group_by(farm) %>% 
  summarise(mean_wt = mean(allE, na.rm = T), sd_wt = sd(allE, na.rm = T))
edgedf %>% 
  group_by(food) %>% 
  summarise(mean_wt = mean(allE, na.rm = T), sd_wt = sd(allE, na.rm = T))
########################################
## 4c. enviro rice comparisons ####
########################################
## enviro rice vs rice named 
my_box_funt(net = riceE, 
            naming = recip.name, 
            my.title  = "enviro rice - complete naming")
wilcox.test(edgedf$riceE[which(edgedf$recip.name == "not.named")], edgedf$riceE[which(edgedf$recip.name == "named")])
wilcox.test(edgedf$riceE[which(edgedf$recip.name == "not.named")], edgedf$riceE[which(edgedf$recip.name == "recip.named")])
wilcox.test(edgedf$riceE[which(edgedf$recip.name == "named")], edgedf$riceE[which(edgedf$recip.name == "recip.named")])

## enviro rice vs freetime 
my_box_funt(net = riceE, 
            naming = freetime, 
            my.title  = "enviro rice - freetime naming")
wilcox.test(edgedf$riceE[which(edgedf$freetime == "not.named")], edgedf$riceE[which(edgedf$freetime == "named")])
wilcox.test(edgedf$riceE[which(edgedf$freetime == "not.named")], edgedf$riceE[which(edgedf$freetime == "recip.named")])
wilcox.test(edgedf$riceE[which(edgedf$freetime == "named")], edgedf$riceE[which(edgedf$freetime == "recip.named")])

## enviro rice vs farm 
my_box_funt(net = riceE, 
            naming = farm, 
            my.title  = "enviro rice - farm naming")
wilcox.test(edgedf$riceE[which(edgedf$farm == "not.named")], edgedf$riceE[which(edgedf$farm == "named")])
wilcox.test(edgedf$riceE[which(edgedf$farm == "not.named")], edgedf$riceE[which(edgedf$farm == "recip.named")])
wilcox.test(edgedf$riceE[which(edgedf$farm == "named")], edgedf$riceE[which(edgedf$farm == "recip.named")])

## enviro rice vs food 
my_box_funt(net = riceE, 
            naming = food, 
            my.title  = "enviro rice - food naming")
wilcox.test(edgedf$riceE[which(edgedf$food == "not.named")], edgedf$riceE[which(edgedf$food == "named")])
wilcox.test(edgedf$riceE[which(edgedf$food == "not.named")], edgedf$riceE[which(edgedf$food == "recip.named")])
wilcox.test(edgedf$riceE[which(edgedf$food == "named")], edgedf$riceE[which(edgedf$food == "recip.named")])

## table of mean and sd weights
edgedf %>% 
  group_by(recip.name) %>% 
  summarise(mean_wt = mean(riceE, na.rm = T), sd_wt = sd(riceE, na.rm = T))
edgedf %>% 
  group_by(freetime) %>% 
  summarise(mean_wt = mean(riceE, na.rm = T), sd_wt = sd(riceE, na.rm = T))
edgedf %>% 
  group_by(farm) %>% 
  summarise(mean_wt = mean(riceE, na.rm = T), sd_wt = sd(riceE, na.rm = T))
edgedf %>% 
  group_by(food) %>% 
  summarise(mean_wt = mean(riceE, na.rm = T), sd_wt = sd(riceE, na.rm = T))
########################################
## 4d. violin plots ####
########################################
labeldf = data.frame(recip.name = rep(c("not.named", "named", "recip.named"), 3), 
                     net = c(rep("a. close contact",3), 
                             rep("b. entire environmental",3), 
                             rep("c. rice environmental",3)), 
                     letter = c("A","B","C",   "A","B","B",   "A","A","B"))

# pdf(file = "plots/violin_plots_cc.pdf", bg = "transparent")
# png(file = "plots/violin_plots_cc.png", bg = "transparent", width = 20, height = 20, units = 'cm', res = 300)
edgedf %>%
  select(id1, id2, recip.name, riceE, allE, cc_weight) %>%
  pivot_longer(cols = c(riceE, allE, cc_weight), names_to = "net", values_to = "edgew") %>%
  mutate(net = case_when(net == "allE" ~ "b. entire environmental",
                         net == "riceE" ~ "c. rice environmental",
                         net == "cc_weight" ~ "a. close contact")) %>%
  group_by(net) %>%
  mutate(edgenorm = edgew/max(edgew, na.rm = T)) %>%
  ggplot(aes(y = edgenorm, x = recip.name)) +
  geom_violin(adjust = 0.5, na.rm = T, trim = T, aes(fill = recip.name, col = recip.name)) +
  geom_boxplot(na.rm = T, width=0.1, color="black", alpha=0.9, outlier.size = 0.5, show.legend = F,
               aes(fill = recip.name)) +
  geom_text(data=labeldf, aes(x=recip.name, y=2, label = letter)) + #for letter labels
  facet_wrap(vars(net))+
  scale_y_continuous(trans = "log10")+
  labs(y = "Edge weight, normalized 0-1 (log 10)") +
  scale_fill_brewer(palette = "Dark2", labels = c("not named", "named", "reciprocated")) +
  scale_color_brewer(palette = "Dark2", labels = c("not named", "named", "reciprocated")) +
  theme(axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(colour = "black")
  )
# dev.off()

########################################
## 5. Look at the outliers ####
########################################
demo <- read.csv("big_edgelist_Mandena1.csv") %>% 
  select(id1, id2, VI_50, gender_match, age_diff, housedist, liveclose)

########################################
## 5a. close contact model  ####
########################################
ol <- edgedf %>% 
  left_join(demo, by = c("id1", "id2")) %>% 
  filter(recip.name == "not.named", !is.na(cc_weight))

# outlier threshold
thres = quantile(ol$cc_weight, .75) +(1.5* IQR(ol$cc_weight))
thres * 20 * 24 * 7 # number of contacts in perfect gps working week
ol <- ol %>% 
  mutate(ccol = ifelse(cc_weight > thres, "outlier", "not"))
table(ol$ccol)

## VI_50
ggplot(ol, aes(x = ccol, y = VI_50)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "VI_50 outliers") 
wilcox.test(ol$VI_50[which(ol$ccol == "not")], ol$VI_50[which(ol$ccol == "outlier")],
            alternative = "less") #x is less than y
ol %>% 
  group_by(ccol) %>% 
  summarise(mean(VI_50), sd(VI_50))

## Gender match
table(ol$gender_match, ol$ccol)
chisq.test(ol$gender_match, ol$ccol)

## liveclose
table(ol$liveclose, ol$ccol)
chisq.test(ol$liveclose, ol$ccol)

## house distance
ggplot(ol, aes(x = ccol, y = housedist)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "housedist outliers") 
wilcox.test(ol$housedist[which(ol$ccol == "not")], ol$housedist[which(ol$ccol == "outlier")],
            alternative = "greater")
ol %>% 
  group_by(ccol) %>% 
  summarise(mean(housedist), sd(housedist))

## check if individuals are 'close' on the naming network just didn't name each other (kinda gwesp)
tmp <- ol %>% filter(ccol == "outlier")

sp = c()
for (i in 1:nrow(tmp)) {
  p1 = tmp$id1[i]
  p2 = tmp$id2[i]
  sd = distances(d_net$n.named, v = p1, to = p2, weights = NA) #calc shortest path on network between the individuals
  sp = c(sp, sd)
}
tmp$sp = sp

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c(">3", "3", "2", "1"))) %>%
  ggplot(aes(x = cc_weight, y = housedist, col = name.distance)) + 
  geom_point() +  
  scale_x_continuous(trans = "log10") +
  theme_minimal() +
  scale_color_brewer(palette = "YlOrRd")

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c("1", "2", "3", ">3"))) %>%
  ggplot(aes(x = name.distance, y = housedist, fill = name.distance)) + 
  geom_boxplot()

########################################
## 5c. enviro all ####
########################################
ol <- edgedf %>% 
  left_join(demo, by = c("id1", "id2")) %>% 
  filter(recip.name == "not.named", !is.na(allE))

# outlier threshold
thres = quantile(ol$allE, .75) +(1.5* IQR(ol$allE))
ol <- ol %>% 
  mutate(aeol = ifelse(allE > thres, "outlier", "not"))
table(ol$aeol)

## VI_50
ggplot(ol, aes(x = aeol, y = VI_50)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "VI_50 outliers") 
wilcox.test(ol$VI_50[which(ol$aeol == "not")], ol$VI_50[which(ol$aeol == "outlier")],
            alternative = "less")
ol %>% 
  group_by(aeol) %>% 
  summarise(mean(VI_50), sd(VI_50))

## Gender match
table(ol$gender_match, ol$aeol)
chisq.test(ol$gender_match, ol$aeol)

## liveclose
table(ol$liveclose, ol$aeol)
chisq.test(ol$liveclose, ol$aeol)

## house distance
ggplot(ol, aes(x = aeol, y = housedist)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "housedist outliers") 
wilcox.test(ol$housedist[which(ol$aeol == "not")], ol$housedist[which(ol$aeol == "outlier")],
            alternative = "greater")
ol %>% 
  group_by(aeol) %>% 
  summarise(mean(housedist), sd(housedist))

## check if individuals are 'close' on the naming network just didn't name each other (kinda gwesp)
tmp <- ol %>% filter(aeol == "outlier")

sp = c()
for (i in 1:nrow(tmp)) {
  p1 = tmp$id1[i]
  p2 = tmp$id2[i]
  sd = distances(d_net$n.named, v = p1, to = p2, weights = NA) #calc shortest path on network between the individuals
  sp = c(sp, sd)
}
tmp$sp = sp

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c(">3", "3", "2", "1"))) %>%
  ggplot(aes(x = allE, y = housedist, col = name.distance)) + 
  geom_point() +  
  scale_x_continuous(trans = "log10") +
  theme_minimal() +
  scale_color_brewer(palette = "YlOrRd")

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c("1", "2", "3", ">3"))) %>%
  ggplot(aes(x = name.distance, y = housedist, fill = name.distance)) + 
  geom_boxplot()

########################################
## 5d. enviro rice ####
########################################
ol <- edgedf %>% 
  left_join(demo, by = c("id1", "id2")) %>% 
  filter(recip.name == "not.named", !is.na(riceE))

# outlier threshold
thres = quantile(ol$riceE, .75) +(1.5* IQR(ol$riceE))
ol <- ol %>% 
  mutate(reol = ifelse(riceE > thres, "outlier", "not"))
table(ol$reol)

## VI_50
ggplot(ol, aes(x = reol, y = VI_50)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "VI_50 outliers") 
wilcox.test(ol$VI_50[which(ol$reol == "not")], ol$VI_50[which(ol$reol == "outlier")],
            alternative = "less")
ol %>% 
  group_by(reol) %>% 
  summarise(mean(VI_50), sd(VI_50))

## Gender match
table(ol$gender_match, ol$reol)
chisq.test(ol$gender_match, ol$reol)

## liveclose
table(ol$liveclose, ol$reol)
chisq.test(ol$liveclose, ol$reol)
fisher.test(ol$liveclose, ol$reol)

## house distance
ggplot(ol, aes(x = reol, y = housedist)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "housedist outliers") 
wilcox.test(ol$housedist[which(ol$reol == "not")], ol$housedist[which(ol$reol == "outlier")],
            alternative = "greater")
ol %>% 
  group_by(reol) %>% 
  summarise(mean(housedist), sd(housedist))

## check if individuals are 'close' on the naming network just didn't name each other (kinda gwesp)
tmp <- ol %>% filter(reol == "outlier")

sp = c()
for (i in 1:nrow(tmp)) {
  p1 = tmp$id1[i]
  p2 = tmp$id2[i]
  sd = distances(d_net$n.named, v = p1, to = p2, weights = NA) #calc shortest path on network between the individuals
  sp = c(sp, sd)
}
tmp$sp = sp

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c(">3", "3", "2", "1"))) %>%
  ggplot(aes(x = riceE, y = housedist, col = name.distance)) + 
  geom_point() +  
  scale_x_continuous(trans = "log10") +
  theme_minimal() +
  scale_color_brewer(palette = "YlOrRd")

tmp %>% 
  mutate(name.distance = case_when(sp == 1 ~ "1",
                                   sp == 2 ~ "2", 
                                   sp == 3 ~ "3",
                                   TRUE ~ ">3")) %>% 
  mutate(name.distance = factor(name.distance, levels = c("1", "2", "3", ">3"))) %>%
  ggplot(aes(x = name.distance, y = housedist, fill = name.distance)) + 
  geom_boxplot()

########################################
## 6. Correlations in centrality ####
########################################
## import centrality tables
name_cent <- read.csv("centrality_name_directed.csv") %>% 
  filter(network == "n.named") %>% select(-network) %>% 
  rename(id = social_netid, name.pagerank = pagerank, name.btwn = btwn, name.strength = totalstrength, name.degree = totaldegree)

cc_cent <- read.csv("centrality_closecontact.csv") %>% 
  rename(cc.degree = med_deg, cc.strength = med_strength, cc.eigen = med_eig, cc.btwn = med_bet)

allE_cent <- read.csv("centrality_enviro.csv") %>% 
  filter(network == "ALL") %>% select(-network) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  rename(id = ID, allE.eigen = EIGEN, allE.btwn = BETWEEN, allE.degree = DEGREE, allE.strength = STRENGTH)

riceE_cent <- read.csv("centrality_enviro.csv") %>% 
  filter(network == "RICE") %>% select(-network) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  rename(id = ID, riceE.eigen = EIGEN, riceE.btwn = BETWEEN, riceE.degree = DEGREE, riceE.strength = STRENGTH)

########################################
## 6a. Top 10% individuals ####
########################################
## naming
nd <- name_cent %>% 
  arrange(desc(name.degree)) %>% 
  pull(id)
ns <- name_cent %>% 
  arrange(desc(name.strength)) %>% 
  pull(id)
nb <- name_cent %>% 
  arrange(desc(name.btwn)) %>% 
  pull(id)
np <- name_cent %>% 
  arrange(desc(name.pagerank)) %>% 
  pull(id)

## close contact model with names
cs <- cc_cent%>% 
  arrange(desc(cc.strength)) %>% 
  pull(id)
cb <- cc_cent %>% 
  arrange(desc(cc.btwn)) %>% 
  pull(id)
ce <- cc_cent %>% 
  arrange(desc(cc.eigen)) %>% 
  pull(id)

## all enviro
as <- allE_cent %>% 
  arrange(desc(allE.strength)) %>% 
  pull(id)
ab <- allE_cent %>% 
  arrange(desc(allE.btwn)) %>% 
  pull(id)
ae <- allE_cent %>% 
  arrange(desc(allE.eigen)) %>% 
  pull(id)

## rice
rs <- riceE_cent %>% 
  arrange(desc(riceE.strength)) %>% 
  pull(id)
rb <- riceE_cent %>% 
  arrange(desc(riceE.btwn)) %>% 
  pull(id)
re <- riceE_cent %>% 
  arrange(desc(riceE.eigen)) %>% 
  pull(id)

top10 <- data.frame(rank = 1:12,
                    name.degree = nd[1:12],
                    name.strength = ns[1:12],
                    name.btwn = nb[1:12],
                    name.page = np[1:12],
                    cc.strength = cs[1:12],
                    cc.btwn = cb[1:12],
                    cc.eigen = ce[1:12],
                    allE.strength = as[1:12],
                    allE.btwn = ab[1:12],
                    allE.eigen = ae[1:12],
                    rice.strength = rs[1:12],
                    rice.btwn = rb[1:12],
                    rice.eigen = re[1:12]
)

## number of networks they are central on regardless of centrality metric
tp10 <- top10 %>% 
  pivot_longer(cols = 2:length(top10), names_to = "metric", values_to = "indiv") %>% 
  separate(metric, c("net", "metric"), sep = "\\.") %>% 
  select(-rank, -metric) %>% 
  count(indiv, net) %>% 
  pivot_wider(names_from = net, values_from = n) %>% 
  mutate(n.networks = 4 - (is.na(rice) + is.na(name) + is.na(cc) + is.na(allE))) %>% 
  arrange(desc(n.networks))
tp10
table(tp10$n.networks)

## by centrality metric
tp10 <- top10 %>%
  pivot_longer(cols = 2:length(top10), names_to = "metric", values_to = "indiv") %>%
  separate(metric, c("net", "metric"), sep = "\\.") %>%
  select(-rank, -net) %>%
  count(indiv, metric) %>%
  pivot_wider(names_from = metric, values_from = n) %>%
  mutate(strength = degree +strength) %>% select(-degree) %>% #tricky thing so naming is degree not strength, but group with strength on other networks
  mutate(n.cent = 3 - (is.na(eigen) + is.na(betweenness) + is.na(strength))) %>%
  arrange(desc(n.cent))
tp10
table(tp10$n.cent)


########################################
## 6b. Spearman's corr of everything  ####
########################################

## heat map 
dfcor <- name_cent %>% 
  left_join(cc_cent, by = "id") %>%
  left_join(allE_cent, by = "id")%>% 
  left_join(riceE_cent, by = "id")%>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  select(-name.strength, -cc.degree, -allE.degree, -riceE.degree) %>% 
  arrange(id) %>% 
  select(-id)

res <- cor(dfcor, method = "spearman")

# write.csv(res, "centrality_correlation_matrix.csv", row.names = F)

col = RColorBrewer::brewer.pal(9, name = "YlOrRd")
heatmap(x = res, col = col,symm = T)

## summary table
diag(res) <- NA
data.frame(min = min(abs(res), na.rm = T),
           max = max(abs(res), na.rm = T), 
           mean = mean(abs(res), na.rm = T),
           sd = sd(abs(res), na.rm = T))

## check for significant correlations
diag(res) <- 1
pres <- data.frame()
for (i in 1:ncol(dfcor)) {
  for (j in 1:ncol(dfcor)) {
    tst = cor.test(dfcor[,i], dfcor[,j], method = "spearman")
    ressig <- data.frame(nd = colnames(dfcor)[i], ck = colnames(dfcor)[j],
                      rho = tst$estimate,
                      pv = tst$p.value)
    pres <- rbind(pres, ressig)
    
  }
}
View(pres)
# write.csv(pres, "centrality_correlations.csv", row.names = F)

pres %>% filter(pv <0.05)
pres <- pres %>% 
  mutate(pvv = case_when(pv <= 0.01 ~ "**",
                         pv <= 0.05 & pv > 0.01 ~ "*",
                         pv > 0.05 & pv <= 0.1 ~ ".",
                         TRUE ~ " ")) %>% 
  select(nd, ck, pvv) %>% 
  pivot_wider(names_from = ck, values_from = pvv)

########################################
## 6c. Spearman's corr of naming with everything else ####
########################################

res1 <- res[c("name.pagerank", "name.btwn", "name.degree"),
            c("cc.strength", "cc.eigen", "cc.btwn" ,
              "allE.eigen", "allE.btwn", "allE.strength",
              "riceE.eigen", "riceE.btwn", "riceE.strength")]
res1
# write.csv(res1, "centrality_correlations_nameVsrest.csv", row.names = F)

#no normalization
heatmap(x = res1, col = col, symm = F, cexRow = 1, cexCol = 0.75, Colv = NA, Rowv = NA, scale = "none")

## row normalization
heatmap(x = res1, col = col, symm = F, cexRow = 1, cexCol = 0.75, Colv = NA, Rowv = NA, scale = "row")

## summary table
data.frame(min = min(abs(res1), na.rm = T),
           max = max(abs(res1), na.rm = T), 
           mean = mean(abs(res1), na.rm = T),
           sd = sd(abs(res1), na.rm = T))

## significant correlations
pres1 <- pres %>% 
  filter(nd %in% c("name.pagerank", "name.btwn", "name.degree")) %>% 
  select(-name.pagerank, -name.btwn, -name.degree)
View(pres1)
# write.csv(pres1, "centrality_correlations_significance_symb.csv", row.names = F)

## Scatter plots of significant correlations (<0.05)  
ggplot(dfcor, aes(x = name.pagerank, y = cc.strength)) + geom_point()
ggplot(dfcor, aes(x = name.pagerank, y = cc.eigen)) + geom_point()
ggplot(dfcor, aes(x = name.pagerank, y = allE.eigen)) + geom_point()

ggplot(dfcor, aes(x = name.btwn, y = cc.btwn)) + geom_point()
ggplot(dfcor, aes(x = name.btwn, y = riceE.eigen)) + geom_point()
ggplot(dfcor, aes(x = name.btwn, y = riceE.btwn)) + geom_point()
ggplot(dfcor, aes(x = name.btwn, y = riceE.strength)) + geom_point()

ggplot(dfcor, aes(x = name.degree, y = cc.strength)) + geom_point()
ggplot(dfcor, aes(x = name.degree, y = cc.btwn)) + geom_point()
ggplot(dfcor, aes(x = name.degree, y = riceE.btwn)) + geom_point()
ggplot(dfcor, aes(x = name.degree, y = riceE.strength)) + geom_point()

## compare metric by metric
## eigen
e <- res[c("name.pagerank", "cc.eigen", "allE.eigen", "riceE.eigen"
), 
c("name.pagerank", "cc.eigen", "allE.eigen", "riceE.eigen")]
diag(e) <- NA
range(e, na.rm = T)
## strength
s <- res[c("name.degree", "cc.strength", "allE.strength", "riceE.strength"), 
         c("name.degree", "cc.strength", "allE.strength", "riceE.strength")]
diag(s) <- NA
range(s, na.rm = T)

## between
b <- res[c("name.btwn", "cc.btwn", "allE.btwn", "riceE.btwn"), 
         c("name.btwn", "cc.btwn", "allE.btwn", "riceE.btwn")]
diag(b) <- NA
range(b, na.rm = T)
