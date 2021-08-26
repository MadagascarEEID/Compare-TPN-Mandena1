## formatting the edges dataframe used in GLM and ERGM (matricies)
library(tidyverse)

########################################
## 1. df of all pairs ####
########################################
## list of everyone that was surveyed
surveyed <- c(paste("A.SNH00", 1:9, sep = ""),paste("A.SNH0", 10:99, sep = ""), paste("A.SNH", 100:176, sep = ""))

## make a non-direction list of all paits
pairs <- data.frame()
for (i in 1:(length(surveyed)-1)) {
  for (j in (i+1):length(surveyed)) {
    tmp <- data.frame(id1 = surveyed[i], id2 = surveyed[j])
    pairs <-rbind(pairs, tmp)
  }
}
pairs$pair <- paste(pairs$id1, pairs$id2, sep = "_")
nrow(pairs); choose(176, 2) #check that is the expected length

## get the list for everyone we have good gps data on
## I am doing this because individuals with less than 240 points might be in GPS derived edge data (prox, vi)
woregps <- read.csv("woregps.csv")[,1]

########################################
## 2. proximity  ####
########################################

### proximity (from GPS data preparation script 8c) 
tmp <- read.csv("prox_3min_Mandena1.csv") %>%
  filter(id1 %in% woregps, id2 %in% woregps, n_simultaneous >= 240)
nrow(tmp);choose(length(woregps),2)
summary(tmp[,3:length(tmp)])

pairs <- pairs %>%
  left_join(tmp, by = c("id1", "id2"))

### coerced proximity (from GPS data preparation script 8d)
tmp <- read.csv("coercedprox_3min_Mandena1.csv") %>%
  filter(id1 %in% woregps, id2 %in% woregps, n_simultaneous >= 240) %>% 
  select(id1, id2, prox) %>% 
  rename(pooled_prox = prox)
nrow(tmp);choose(length(woregps),2)
summary(tmp[,3:length(tmp)])

###
pairs <- pairs %>%
  left_join(tmp, by = c("id1", "id2"))

########################################
## 3. VI  ####
########################################
#### VI_95 and VI_50 
tmp <- read.csv("vi_all_Mandena1.csv") %>% 
  filter(id1 %in% woregps, id2 %in% woregps)  
nrow(tmp);choose(length(woregps),2)
summary(tmp[,3:length(tmp)])

pairs <- pairs %>%
  left_join(tmp, by = c("id1", "id2"))

#### VI day and night
tmp <- read.csv("vi_nightday_Mandena1.csv") %>% 
  filter(id1 %in% woregps, id2 %in% woregps)
nrow(tmp);choose(length(woregps),2)
summary(tmp[,3:length(tmp)]) #vi 50 night max is 50.19 might want to set to 50

pairs <- pairs %>%
  left_join(tmp, by = c("id1", "id2"))

### from the simultaneous wear VI's
tmp <- read.csv("vi_cowear_Mandena1.csv") %>% 
  rename(VI_50_cowear = VI_50, VI_95_cowear = VI_95) %>% 
  filter(id1 %in% woregps, id2 %in% woregps)
nrow(tmp);choose(length(woregps),2)
summary(tmp[,3:length(tmp)])

pairs <- pairs %>% 
  left_join(tmp, by = c("id1", "id2"))

########################################
## 4. naming  ####
########################################

## count number of times named the same person, including for which question
## and mark mutual edges (recip.naming)
library(igraph)
nd_net <- readRDS("network_name_nondirected.RDS")
d_net <- readRDS("network_name_directed.RDS")

## n.named and recip.named
# n.named from undirected network (edges summed)
q <- as_data_frame(nd_net$n.named, "edges") %>% 
  select(from, to, weight) %>% 
  rename(id1 = from, id2 = to, n.named = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## q1
q <- as_data_frame(nd_net$freetime, "edges") %>% 
  rename(id1 = from, id2 = to, freetime = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## q2
q <- as_data_frame(nd_net$farm_helpyou, "edges") %>% 
  rename(id1 = from, id2 = to, farm_helpyou = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## q3
q <- as_data_frame(nd_net$farm_youhelp, "edges") %>% 
  rename(id1 = from, id2 = to, farm_youhelp = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## q4
q <- as_data_frame(nd_net$food_helpyou, "edges") %>% 
  rename(id1 = from, id2 = to, food_helpyou = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## q5
q <- as_data_frame(nd_net$food_youhelp, "edges") %>% 
  rename(id1 = from, id2 = to, food_youhelp = weight)

pairs <- pairs %>%  
  left_join(q, by = c("id1", "id2"))

## recip.named from mutual network
E(d_net$n.named)$recip.name = which_mutual(d_net$n.named)
recip <- as_data_frame(d_net$n.named) %>% 
  select(from, to, recip.name) %>% 
  mutate(recip.name = ifelse(recip.name, 1, 0))

recip <- recip %>%
  select(from, to, recip.name)  

pairs <- recip %>% 
  rename(id1 = from, id2 = to) %>%  
  right_join(pairs, by = c("id1", "id2")) %>% 
  mutate(recip.name = ifelse(n.named > 0 & is.na(recip.name), 0, recip.name),
         recip.name = case_when(recip.name == 0 ~ "named",
                                recip.name == 1 ~ "recip.named",
                                TRUE ~ "not.named"))

detach(package:igraph)

pairs %>% 
  mutate(ck = freetime + farm_helpyou + farm_youhelp + food_helpyou + food_youhelp) %>% 
  filter(n.named != ck)
table(pairs$recip.name)
rm(q, nd_net, d_net)

## make the NA on the questions into 0s (we know they didn't name eachother...)
pairs <- pairs %>% 
  mutate(across(c(n.named, freetime, farm_helpyou, farm_youhelp, food_helpyou, food_youhelp), 
                ~replace_na(.x, 0)))

########################################
## 5. demographic data  ####
########################################

demo <- read.csv("coded_demodat.csv") %>% 
  select(social_netid, gender, age)

demo <- demo %>% 
  mutate(age = ifelse(is.na(age), round(mean(demo$age, na.rm = T)), age))

demo.match <- pairs %>% 
  select(id1, id2) %>%
  left_join(demo, by = c("id1" = "social_netid")) %>% 
  left_join(demo, by = c("id2" = "social_netid")) %>% 
  mutate(gender_match = ifelse(gender.x == gender.y, "same", "different")) %>% # gender match
  mutate(age_diff = abs(age.x-age.y))  %>% # age difference
  select(id1, id2, gender_match, age_diff)

pairs <- pairs %>% 
  left_join(demo.match, by = c("id1", "id2"))
rm(demo.match)

########################################
## 6. house distance and live close  ####
########################################

### house distance was calculated using this code. hc = df of participants house location (utm)
# house.dist <- data.frame()
# for (i in 1:(nrow(hc)-1)) {
#   for (j in 2:nrow(hc)) { #typo here it should be (i+1):nrow(hc), results in making more comparisons than needed
#     send <- hc[i,]
#     reci <- hc[j,] 
#     out <- data.frame(sender = send$social_netid, 
#                       receiver = reci$social_netid, 
#                       dist = sqrt((send$xcoord - reci$xcoord)^2 + (send$ycoord - reci$ycoord)^2))
#     house.dist <- rbind(house.dist, out)
#     
#     
#   }
# }
tmp <- read.csv("housedist_Mandena12.csv") %>% 
  select(-X) %>% 
  rename(id1 = sender, id2 = receiver, housedist = dist) %>% 
  mutate(id1 = str_replace(id1, "-", ".")) %>% 
  mutate(id2 = str_replace(id2, "-", ".")) 

pairs<- pairs %>% 
  left_join(tmp, by = c("id1", "id2")) %>% 
  mutate(liveclose = as.factor(ifelse(housedist >25, 1, 0))) #this was 0, 1 before I think that was backwards though!
summary(pairs$housedist)
summary(pairs$liveclose)
rm(tmp, surveyed)

########################################
## 7. save df ####
########################################

write.csv(pairs, "big_edgelist_Mandena1.csv", row.names = F)

########################################
## 8. make matricies ####
########################################
require(igraph)

## get just the gps wearer edges
eedges <- pairs %>% 
  filter(id1 %in% woregps, id2 %in% woregps) %>% 
  select(id1, id2, prox, pooled_prox, starts_with("VI"), housedist, liveclose, n.named, recip.name) %>% 
  mutate(liveclose = as.numeric(liveclose),
         recip.name = ifelse(recip.name != "recip.named", 0, 1))

## build a matrix for each variable
g <- graph_from_data_frame(eedges, vertices = woregps, directed = FALSE)

colnames(eedges)
m.prox = as_adjacency_matrix(g, attr = "prox") %>% as.matrix()
m.pprox = as_adjacency_matrix(g, attr = "pooled_prox") %>% as.matrix()
m.vi95 = as_adjacency_matrix(g, attr = "VI_95") %>% as.matrix()
m.vi50 = as_adjacency_matrix(g, attr = "VI_50") %>% as.matrix()
m.vi95n = as_adjacency_matrix(g, attr = "VI_95_night") %>% as.matrix()
m.vi50n = as_adjacency_matrix(g, attr = "VI_50_night") %>% as.matrix()
m.vi95d = as_adjacency_matrix(g, attr = "VI_95_day") %>% as.matrix()
m.vi50d = as_adjacency_matrix(g, attr = "VI_50_day") %>% as.matrix()
m.vi95c = as_adjacency_matrix(g, attr = "VI_95_cowear") %>% as.matrix()
m.vi50c = as_adjacency_matrix(g, attr = "VI_50_cowear") %>% as.matrix()
m.hdist = as_adjacency_matrix(g, attr = "housedist") %>% as.matrix()
m.lc = as_adjacency_matrix(g, attr = "liveclose") %>% as.matrix()
m.nnamed = as_adjacency_matrix(g, attr = "n.named") %>% as.matrix()
m.recip.named = as_adjacency_matrix(g, attr = "recip.name") %>% as.matrix()


## make a list of those matrices to save
ergmlist <- list(m.prox, 
                 m.pprox,
                 m.vi95, 
                 m.vi50, 
                 m.vi95n, 
                 m.vi50n, 
                 m.vi95d, 
                 m.vi50d, 
                 m.vi95c,
                 m.vi50c,
                 m.hdist, 
                 m.lc,
                 m.nnamed,
                 m.recip.named)

names(ergmlist) = c("observed_prox", 
                    "pooled_prox", 
                    "VI_95", 
                    "VI_50", 
                    "VI_95_night", 
                    "VI_50_night", 
                    "VI_95_day",
                    "VI_50_day",
                    "VI_95_cowear",
                    "VI_50_cowear",
                    "housedist",
                    "liveclose",
                    "n.named",
                    "recip.named")

saveRDS(ergmlist, "ergmmatrices_Mandena1.RDS")
