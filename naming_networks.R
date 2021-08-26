## Naming network
require(tidyverse)
require(igraph)

########################################
## 1. Load the data  ####
########################################

## the edge df interviewee = sender, person they named = receiver
### network_question_number = 1-5, see table 1 in ms
### this data was collected during the survey, everyone named was later assigned an uid
edgedf <- read.csv(file = "NameTable_Mandena12.csv") %>% 
  select(social_netid, named_person_social_netid, network_question_number) %>% 
  rename(sender = social_netid, receiver = named_person_social_netid) %>% 
  mutate(across(c("sender", "receiver"), ~str_replace(.x, "-", "."))) # in the social_netid's change the "-" to a "."

## list of indivs with good gps data
woregps <- read.csv(file = "woregps.csv")[,1]

## filter edge df to just gps wearing senders/recievers
gpsedges <- edgedf %>% 
  filter(sender %in% woregps, receiver %in% woregps)

########################################
## 2. Make question specific networks ####
########################################

## check no errors in data resulting in naming a person more than once per questions
gpsedges %>% count(network_question_number, sender, receiver) %>% filter(n >1)

## count number of times named the same person, including for which question
d_names <- gpsedges %>% 
  mutate(did_name = 1) %>% 
  pivot_wider(id_col = c("sender", "receiver"),
              names_from = network_question_number, 
              names_prefix = "q", 
              values_from = did_name,
              values_fill = 0) %>%
  mutate(n.named = q1+q2+q3+q4+q5, 
         n.helping = q2+q3+q4+q5, #never describe this network. check that and remove
         q23 = q2+q3, 
         q45 = q4 + q5, 
         q24 = q2+q4, #remove
         q35 = q3+q5) #remove
summary(d_names)

## function to make each network  
make_networks <- function(x, question, directed, recip.edges){
  x <- x %>% 
    rename("question" = question) %>% 
    select(sender, receiver, question) %>% 
    filter(question >0) #to remove edges that werent observed
  
  
  g <- graph_from_data_frame(x, directed = directed, vertices = woregps) #verticies must equal all gps wearers b/c disconnected graph
  E(g)$weight <- E(g)$question #weight by number of times named
  if(recip.edges){ #detect reciprocal naming
    E(g)$mutual <- ifelse(which_mutual(g), "blue", "grey50")
  }
  return(g)
}

## make each DIRECTED network
d_all <- make_networks(d_names, "n.named", directed = TRUE, recip.edges = TRUE)
d_help <- make_networks(d_names, "n.helping", directed = TRUE, recip.edges = TRUE) #remove?
d_q1 <- make_networks(d_names, "q1", directed = TRUE, recip.edges = TRUE)
d_q2 <- make_networks(d_names, "q2", directed = TRUE, recip.edges = TRUE)
d_q3 <- make_networks(d_names, "q3", directed = TRUE, recip.edges = TRUE)
d_q4 <- make_networks(d_names, "q4", directed = TRUE, recip.edges = TRUE)
d_q5 <- make_networks(d_names, "q5", directed = TRUE, recip.edges = TRUE)
d_farm <- make_networks(d_names, "q23", directed = TRUE, recip.edges = TRUE)
d_food <- make_networks(d_names, "q45", directed = TRUE, recip.edges = TRUE)
d_helpyou <- make_networks(d_names, "q24", directed = TRUE, recip.edges = TRUE) #remove
d_youhelp <- make_networks(d_names, "q35", directed = TRUE, recip.edges = TRUE) #remove

## make those into UNDIRECTED networks
nd_all <- as.undirected(d_all, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_help <- as.undirected(d_help, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore")) #REMOVE?
nd_q1 <- as.undirected(d_q1, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_q2 <- as.undirected(d_q2, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_q3 <- as.undirected(d_q3, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_q4 <- as.undirected(d_q4, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_q5 <- as.undirected(d_q5, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_farm <- as.undirected(d_farm, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_food <- as.undirected(d_food, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore"))
nd_helpyou <- as.undirected(d_helpyou, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore")) #REMOVE
nd_youhelp <- as.undirected(d_youhelp, mode = "collapse", edge.attr.comb = list(weight="sum", "ignore")) #REMOVE

## save those as lists

saveRDS(list(n.named = d_all, 
             help = d_help, 
             freetime = d_q1, 
             farm_helpyou = d_q2, 
             farm_youhelp = d_q3,
             food_helpyou = d_q4, 
             fodd_youhelp = d_q5, 
             farm = d_farm, 
             food = d_food, 
             helpyou = d_helpyou, 
             youhelp = d_youhelp), file = "network_name_directed.RDS") ## remove help
saveRDS(list(n.named = nd_all, 
             help = nd_help, 
             freetime = nd_q1, 
             farm_helpyou = nd_q2, 
             farm_youhelp = nd_q3,
             food_helpyou = nd_q4, 
             food_youhelp = nd_q5, 
             farm = nd_farm, 
             food = nd_food, 
             helpyou = nd_helpyou, 
             youhelp = nd_youhelp), file = "network_name_nondirected.RDS") #remove help
########################################
## 3. network wide summaries ####
########################################

## make a function for calculating the wide stats 
wide_function <- function(x, directed){
  connected_strong = is.connected(x, mode = "strong")
  connected_weak = is.connected(x, mode = "weak")
  diameter = diameter(x)
  avg_distance = average.path.length(x)
  
  density = graph.density(x)
  transitivity_global = transitivity(x, type = "global")
  if (directed) { 
    reciprocity = reciprocity(x) #CAN calculate reciprocity on directed network
    louvain_modularity = NA #CANNOT calculate modularity on directed network
  } else{
    reciprocity = NA #CANNOT calculate reciprocity on directed network
    louvain_modularity = modularity(cluster_louvain(x)) #CAN calculate modularity on directed network
  }
  return(data.frame(connected_strong, 
                    connected_weak, 
                    diameter, 
                    avg_distance, 
                    reciprocity, 
                    density, 
                    transitivity_global, 
                    louvain_modularity)
  )
  
}

## directed
widestats.d <- rbind(wide_function(d_all, directed = TRUE),
                     wide_function(d_q1, directed = TRUE),
                     wide_function(d_q2, directed = TRUE),
                     wide_function(d_q3, directed = TRUE),
                     wide_function(d_q4, directed = TRUE),
                     wide_function(d_q5, directed = TRUE),
                     wide_function(d_farm, directed = TRUE),
                     wide_function(d_food, directed = TRUE)) 
widestats.d <- widestats.d %>% 
  mutate(network = c("n.named", 
                     "free.time", 
                     "help.you.farm", 
                     "you.help.farm", 
                     "help.you.food", 
                     "you.help.food", 
                     "farm",
                     "food")) %>%
  select(network, everything())
widestats.d

## undirected
widestats.nd <- rbind(wide_function(nd_all, directed = FALSE),
                      wide_function(nd_q1, directed = FALSE),
                      wide_function(nd_q2, directed = FALSE),
                      wide_function(nd_q3, directed = FALSE),
                      wide_function(nd_q4, directed = FALSE),
                      wide_function(nd_q5, directed = FALSE),
                      wide_function(nd_farm, directed = FALSE),
                      wide_function(nd_food, directed = FALSE))
widestats.nd <- widestats.nd %>% 
  mutate(network = c("n.named",
                     "free.time", 
                     "help.you.farm", 
                     "you.help.farm", 
                     "help.you.food", 
                     "you.help.food", 
                     "farm",
                     "food")) %>% 
  select(network, everything())
widestats.nd

## save these
write.csv(widestats.d, file = "widestats_name_directed.csv", row.names = F)
write.csv(widestats.nd, file = "widestats_name_nondirected.csv", row.names = F)

# Look at the size of the components:    
groups(components(d_all))
groups(components(d_q1))
groups(components(d_q2)) 
groups(components(d_q3))
groups(components(d_q4))
groups(components(d_q5))
groups(components(d_farm))
groups(components(d_food))

########################################
## 4. Centrality ####
########################################

cent_function <- function(x, directed){
  ## total degree
  cent_tab <- data.frame(degree(x, mode = "total")) %>% 
    rownames_to_column("social_netid") %>% 
    rename(totaldegree = starts_with("degree"))
  
  ## total strength
  cent_tab <- data.frame(strength(x, mode = "total", weights = E(x)$weight)) %>% 
    rownames_to_column("social_netid") %>% 
    rename(totalstrength = starts_with("strength")) %>% 
    right_join(cent_tab)
  
  ## betweenness
  cent_tab <- data.frame(betweenness(x, normalized = TRUE, weights = 1-E(x)$weight/(max(E(x)$weight)+1))) %>% 
    rownames_to_column("social_netid") %>% 
    rename(btwn = starts_with("between")) %>% 
    right_join(cent_tab)
  
  if (directed) { ## check if remove in/out deg or strength REMOVE
    # ## in degree
    # cent_tab <- data.frame(degree(x, mode = "in")) %>% 
    #   rownames_to_column("social_netid") %>% 
    #   rename(indegree = starts_with("degree")) %>% 
    #   right_join(cent_tab)
    # ## out degree
    # cent_tab <- data.frame(degree(x, mode = "out")) %>% 
    #   rownames_to_column("social_netid") %>% 
    #   rename(outdegree = starts_with("degree")) %>% 
    #   right_join(cent_tab)
    # ## in strength
    # cent_tab <- data.frame(strength(x, mode = "in", weights = E(x)$weight)) %>% 
    #   rownames_to_column("social_netid") %>% 
    #   rename(instrength = starts_with("strength")) %>% 
    #   right_join(cent_tab)
    # ## out strength
    # cent_tab <- data.frame(strength(x, mode = "out", weights = E(x)$weight)) %>% 
    #   rownames_to_column("social_netid") %>% 
    #   rename(outstrength = starts_with("strength")) %>% 
    #   right_join(cent_tab)
    ## pagerank
    cent_tab <- data.frame(page_rank(x, algo = "arpack", directed = TRUE, weights = E(x)$weight)$vector) %>% 
      rownames_to_column("social_netid") %>% 
      rename(pagerank = starts_with("page")) %>%
      mutate(pagerank = scale(pagerank) + (-min(scale(pagerank)))) %>% #scale pagerank and put on a positive scale
      right_join(cent_tab)
    
  } else{
    ## eigenvector
    cent_tab <- data.frame(eigen_centrality(x, directed = FALSE, weights = NA)$vector) %>% 
      rownames_to_column("social_netid") %>% 
      rename(eigen = starts_with("eigen")) %>% 
      right_join(cent_tab)
    ## make undirected data.frame
    
  }
  return(cent_tab)
}

## directed
cent.d <- rbind(cent_function(d_all, directed = TRUE),
                cent_function(d_q1, directed = TRUE),
                cent_function(d_q2, directed = TRUE),
                cent_function(d_q3, directed = TRUE),
                cent_function(d_q4, directed = TRUE),
                cent_function(d_q5, directed = TRUE),
                cent_function(d_farm, directed = TRUE),
                cent_function(d_food, directed = TRUE))
cent.d <- cent.d %>% 
  mutate(network = c(rep("n.named", 123), 
                     rep("free.time", 123), 
                     rep("help.you.farm", 123), 
                     rep("you.help.farm", 123), 
                     rep("help.you.food", 123), 
                     rep("you.help.food", 123), 
                     rep("farm", 123),
                     rep("food", 123))) %>% 
  select(network, everything())

## undirected
cent.nd <- rbind(cent_function(nd_all, directed = FALSE),
                 cent_function(nd_q1, directed = FALSE),
                 cent_function(nd_q2, directed = FALSE),
                 cent_function(nd_q3, directed = FALSE),
                 cent_function(nd_q4, directed = FALSE),
                 cent_function(nd_q5, directed = FALSE),
                 cent_function(nd_farm, directed = FALSE),
                 cent_function(nd_food, directed = FALSE))
cent.nd <- cent.nd %>% 
  mutate(network = c(rep("n.named", 123), 
                     rep("free.time", 123), 
                     rep("help.you.farm", 123), 
                     rep("you.help.farm", 123), 
                     rep("help.you.food", 123), 
                     rep("you.help.food", 123), 
                     rep("farm", 123),
                     rep("food", 123))) %>% 
  select(network, everything())

## save the scores 
write.csv(cent.d, file = "centrality_name_directed.csv", row.names = F)
write.csv(cent.nd, file = "centrality_name_nondirected.csv", row.names = F)
