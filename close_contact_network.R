## R version 4.0.5 (2021-03-31)

## impute close contact edges using an
## ERGM to determine if there is a close contact edge on the imputed network
## and a GLM to assign edge weights on the network
################################################################################

########################################
## 1. Data preparation ####
########################################
## see close_contact_predictors 

woregps <- read.csv("woregps.csv")[,1]

########################################
## 2. ERGM ~ edge probability ####
########################################
library(tidyverse)
library(network)
library(ergm)

m <- readRDS("ergmmatrices_Mandena1.RDS")
names(m)

dem_data <- read.csv("coded_demodat.csv") %>% 
  select(social_netid, age, gender) %>% 
  filter(social_netid %in% woregps) %>% 
  mutate(age = ifelse(is.na(age), mean(age, na.rm = T), age))

########################################
## 2.1 create binary edges ####
########################################
## In order to run an ERGM for imputation of missing edges, the edges must be binary. 
### Therefore, we must choose a threshold for weighted prox and create a new column in which 
### all values higher than that threshold are 1. This must be done twice-- once for pooled prox 
### (used to fit the ERGM) and once for observed prox (used to simulate networks from the ERGM). 
### Pooled prox and observed prox are highly correlated (see plots below) but do not have an exact 1:1 
### relationship-- therefore, we will use a percentile to determine the threshold value.

## observed prox (stored as a matrix)
t1 <- median(m$observed_prox[which(m$observed_prox > 0)])
observed_prox_binary_matrix <- ifelse(m$observed_prox > t1, 1, 0)
### number of missing edges
length(which(is.na(observed_prox_binary_matrix)))
### remove edge with self
diag(observed_prox_binary_matrix) <- NA
### then calculate network density 
length(which(observed_prox_binary_matrix == 1))/(length(which(observed_prox_binary_matrix == 1))+length(which(observed_prox_binary_matrix == 0)))
## threshold in contacts per day
t1*20*24

## pooled prox (stored as a matrix)
t2 <- median(m$pooled_prox[which(m$pooled_prox > 0)])
pooled_prox_binary_matrix <- ifelse(m$pooled_prox > t2, 1, 0)
### number of missing edges
length(which(is.na(pooled_prox_binary_matrix)))
### remove edge with self
diag(pooled_prox_binary_matrix) <- NA
### then calculate network density 
length(which(pooled_prox_binary_matrix == 1))/(length(which(pooled_prox_binary_matrix == 1))+length(which(pooled_prox_binary_matrix == 0)))

## threshold in contacts per day
t2*20*24
########################################
## 2.2 make both networks ####
########################################
## observed prox network 
observed_prox_new_binary_net = as.network(observed_prox_binary_matrix , directed=FALSE , loops = FALSE , matrix.type = "a" , ignore.eval=FALSE, names.eval="prox_binary")
#add attributes
observed_prox_new_binary_net%v%"gender" <- dem_data$gender
observed_prox_new_binary_net%v%"age" <- dem_data$age

## pooled prox network
pooled_prox_new_binary_net = as.network(pooled_prox_binary_matrix , directed=FALSE , loops = FALSE , matrix.type = "a" , ignore.eval=FALSE, names.eval="prox_binary")
#add attributes
pooled_prox_new_binary_net%v%"gender" <- dem_data$gender
pooled_prox_new_binary_net%v%"age" <- dem_data$age

########################################
## 2.3 test/train datasets ####
########################################
## split observed prox in to training (80%) and testing (20%)
set.seed(246)
diag(observed_prox_binary_matrix) <- 0 # add edge with self back in
full.dat <- reshape2::melt(observed_prox_binary_matrix) %>% 
  rename(id1 = Var1, id2 = Var2, prox_binary = value) %>% 
  drop_na()
train.index <- sample(1:nrow(full.dat), 
                      0.8 * nrow(full.dat))
train.dat <- full.dat[train.index,]
test.dat <- full.dat[-train.index,]
length(which(test.dat$prox_binary == 1))
length(which(test.dat$prox_binary == 0))
length(which(is.na(test.dat$prox_binary)))

train.dat.matrix<- reshape2::dcast(train.dat, id1 ~ id2, value.var = "prox_binary")
row.names(train.dat.matrix) <-  train.dat.matrix$id1
train.dat.matrix <- train.dat.matrix %>% 
  select(-id1) %>% 
  as.matrix

#create the network
train.dat.net = as.network(train.dat.matrix , directed=FALSE , loops = FALSE , matrix.type = "a" , 
                           ignore.eval=FALSE, names.eval="prox_binary")
#add attributes
train.dat.net%v%"gender" <- dem_data$gender
train.dat.net%v%"age" <- dem_data$age

########################################
## 2.4 sens/spec/accuracy function ####
########################################
accuracy <- function(sim_network) {
  
  accuracy <- list()
  specificity <- list()
  sensitivity <- list()
  
  set.seed(246)
  
  for(i in 1:length(sim_network)) {
    #make each matrix into a df
    output <- as.matrix(sim_network[[i]]) 
    edgelist <- reshape2::melt(output)
    colnames(edgelist) <- c("id1", "id2", "predicted")
    
    #join df to test.dat
    accuracy_check <- left_join(test.dat, edgelist, by = c("id1", "id2")) 
    
    #test if prediction is accurate
    accuracy_check <- accuracy_check %>%
      mutate(accurate_prediction = ifelse(prox_binary == predicted, "Yes", "No"))
    
    contab <- table(accuracy_check$predicted, accuracy_check$prox_binary)
    sensitivity[[i]] <- contab[2,2] / (contab[2,2] + contab[2,1])
    specificity[[i]] <- contab[1,1] / (contab[1,1] + contab[1,2])
    accuracy[[i]] <- (contab[2, 2] + contab[1,1])/sum(contab)
  }
  
  overall_df <- as.data.frame(do.call(cbind, accuracy))
  overall_df2 <- as.data.frame(t(overall_df))
  
  specificity_df <- as.data.frame(do.call(cbind, specificity))
  specificity_df2 <- as.data.frame(t(specificity_df))
  
  sensitivity_df <- as.data.frame(do.call(cbind, sensitivity))
  sensitivity_df2 <- as.data.frame(t(sensitivity_df))
  
  return(list("accuracy" = mean(overall_df2$V1), 
              "specificity" = mean(specificity_df2$V1), 
              "sensitivity" = mean(sensitivity_df2$V1)))
}

########################################
## 2.5 proximity plots ####
########################################
## observed prox all values
reshape2::melt(m$observed_prox, na.rm = T) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(binwidth = 0.005) + 
  xlab("Observed Prox")

## observed prox non-zero values
reshape2::melt(m$observed_prox, na.rm = T) %>% 
  filter(value >0) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(binwidth = 0.005) + 
  xlab("Nonzero Observed Prox") + 
  geom_vline(xintercept=t1, linetype = "dashed") #binary cutoff line

## pooled prox all values
reshape2::melt(m$pooled_prox) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(binwidth = 0.005) + 
  xlab("Pooled Prox")

## pooled prox non-zero values
reshape2::melt(m$pooled_prox) %>% 
  filter(value >0) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(binwidth = 0.005) + 
  xlab("Nonzero Pooled Prox") + 
  geom_vline(xintercept=t2, linetype = "dashed") #binary cutoff line

## observed vs pooled prox
reshape2::melt(m$observed_prox, na.rm = T) %>% 
  left_join(reshape2::melt(m$pooled_prox), by = c("Var1", "Var2")) %>% 
  ggplot(aes(x = value.x, y = value.y)) + 
  geom_point() +
  xlab("Observed Prox") + ylab("Pooled prox") + 
  geom_abline() +
  geom_hline(yintercept = t2, linetype = "dashed") + #pooled prox binary cutoff
  geom_vline(xintercept = t1, linetype = "dashed") #observed prox binary cutoff

reshape2::melt(m$observed_prox, na.rm = T) %>% 
  left_join(reshape2::melt(m$pooled_prox), by = c("Var1", "Var2")) %>% 
  lm(value.x ~ value.y, data = .)

########################################
## 2.6 ergm function ####
########################################
set.seed(246)
myergm <- function(mymodel, modname){
  mod <- ergm(mymodel, 
              control = control.ergm(MCMLE.maxit = 10000 , 
                                     MCMC.burnin = 70000, MCMC.interval = 1000, 
                                     parallel=6, parallel.type = "PSOCK"))
  print(summary(mod))
  
  ## assess goodness of fit (increases computational time a lot!)
  # mod.gof <- gof(mod, GOF = ~model)
  # print(mod.gof)
  # plot(mod.gof, main = modname)
  
  saveRDS(mod, file = paste0(modname, ".RDS")) #if want to save each model
  
  sim_test <- simulate(mod,
                       nsim = 500, seed = NULL,
                       burnin = 1000, interval = 1000,
                       basis = train.dat.net,
                       verbose = FALSE)
  ssa = accuracy(sim_test)
  
  return(list(model = mod,
              accuracy = ssa$accuracy,
              specificity = ssa$specificity,
              sensitivity = ssa$sensitivity))
}

########################################
## 2.7 ergm selection ####
########################################
## general plan is to figure out which vi method, 
### then other location info vs naming, then age and gender, then structural terms
### and since we plan to compare cc network to naming networks later on, models without 
### naming will also be considered

## edges model ####
##	pooled_prox ~ edges
mymod <- myergm((pooled_prox_new_binary_net ~ edges), "ergm_models/m0")
ergm_selection <- data.frame(mod = "m0", aic = summary(mymod$model)$aic, 
                             accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)

## VI models ####
## pooled_prox ~ edges + vi_95
mymod <- myergm((pooled_prox_new_binary_net ~ edges + edgecov(m$VI_95, attrname = "vi_95")), 
                "ergm_models/m1")
tmp <- data.frame(mod = "m1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + vi_50
mymod <- myergm((pooled_prox_new_binary_net ~ edges + edgecov(m$VI_50, attrname = "vi_50")), 
                "ergm_models/m2")
tmp <- data.frame(mod = "m2", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night")),
                "ergm_models/m3")
tmp <- data.frame(mod = "m3", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## night/day best, loc vs naming ####
## pooled_prox ~ edges + day_vi95 + night_vi95 + house distance
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist, attrname = "housedist")),
                "ergm_models/m4")
tmp <- data.frame(mod = "m4", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + house distance*liveclose
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist*m$liveclose, attrname = "house_close")),
                "ergm_models/m5")
tmp <- data.frame(mod = "m5", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named")),
                "ergm_models/m6")
tmp <- data.frame(mod = "m6", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + recip.named
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$recip.named, attrname = "recip.named")),
                "ergm_models/m7")
tmp <- data.frame(mod = "m7", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named")),
                "ergm_models/m8")
tmp <- data.frame(mod = "m8", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## combo of those vs adding age and gender ####
## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named + housedist
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   edgecov(m$housedist, attrname = "housedist")),
                "ergm_models/m9")
tmp <- data.frame(mod = "m9", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named + gender
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   nodematch("gender")),
                "ergm_models/m10")
tmp <- data.frame(mod = "m10", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named + age difference
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   absdiff("age" , pow=1)),
                "ergm_models/m11")
tmp <- data.frame(mod = "m11", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## naming + housedist are fairly similar <1 aic so just use naming 
## gwesp fam ####
## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named + gwesp
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   gwesp(0, fixed = TRUE)),
                "ergm_models/m12.0")
tmp <- data.frame(mod = "m12.0", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   gwesp(1, fixed = TRUE)),
                "ergm_models/m12.1")
tmp <- data.frame(mod = "m12.1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   gwesp(0.5, fixed = TRUE)),
                "ergm_models/m12.05")
tmp <- data.frame(mod = "m12.05", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

##check gof 
gof12.1 <- gof(readRDS("ergm_models/m12.1.RDS"))
gof12.1 #all parameters poor
plot(gof12.1) #esp plot is skewed try adding gwnsp
saveRDS(gof12.1, "ergm_models/gof12.1.RDS")
rm(gof12.1)

## pooled_prox ~ edges + day_vi95 + night_vi95 + n.named + recip.named + gwesp + gwnsp
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m13.1")
tmp <- data.frame(mod = "m13.1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## final comparisons ####
## try last model with house dist
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   edgecov(m$housedist, attrname = "housedist") +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m14.1")
tmp <- data.frame(mod = "m14.1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## check gof 
gof14.1 <- gof(readRDS("ergm_models/m14.1"))
gof14.1
plot(gof14.1)
saveRDS(gof14.1, "ergm_models/gof14.1.RDS")

## put everything in the model, compare full vi95 and night/day
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   edgecov(m$housedist*m$liveclose, attrname = "houseclose") + 
                   nodematch("gender") + absdiff("age" , pow=1) +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m15.1")
tmp <- data.frame(mod = "m15.1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## vi_95 instead
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95, attrname = "vi_95") +
                   edgecov(m$n.named, attrname = "n.named") + edgecov(m$recip.named, attrname = "recip.named") +
                   edgecov(m$housedist, attrname = "housedist") +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m16.1")
tmp <- data.frame(mod = "m16.1", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## models without naming infor (n.named and recip.named) ####
### from m4 and m5 know housedist is housedist*liveclose are ~ the same so compare both with demographics
## pooled prox ~ vi 95 day + vi 95 night + housedist*liveclose + gender match + age diff
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist*m$liveclose, attrname = "houseclose") + 
                   nodematch("gender") + absdiff("age" , pow=1) ),
                "ergm_models/m17")
tmp <- data.frame(mod = "m17", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled prox ~ vi 95 day + vi 95 night + housedist + gender match + age diff
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist, attrname = "housedist") + 
                   nodematch("gender") + absdiff("age" , pow=1) ),
                "ergm_models/m18")
tmp <- data.frame(mod = "m18", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## still about the same, add in gwesp
## pooled prox ~ vi 95 day + vi 95 night + housedist*liveclose + gender match + age diff +gwesp(1, fixed)
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist*m$liveclose, attrname = "houseclose") + 
                   nodematch("gender") + absdiff("age" , pow=1) +
                   gwesp(1, fixed = TRUE)),
                "ergm_models/m19")
tmp <- data.frame(mod = "m19", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled prox ~ vi 95 day + vi 95 night + housedist + gender match + age diff +gwesp(1, fixed)
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist, attrname = "housedist") + 
                   nodematch("gender") + absdiff("age" , pow=1) +
                   gwesp(1, fixed = TRUE)),
                "ergm_models/m20")
tmp <- data.frame(mod = "m20", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## They are still about the same, check fit then add gwnsp
##check gof model 19
gof19 <- gof(readRDS("ergm_models/m19.RDS"))
gof19 #all parameters poor
plot(gof19) #esp plot is skewed try adding gwnsp
saveRDS(gof19, "ergm_models/gof19.RDS")
rm(gof19)

## model 20
gof20 <- gof(readRDS("ergm_models/m20.RDS"))
gof20 #all parameters poor
plot(gof20) #esp plot is skewed try adding gwnsp
saveRDS(gof20, "ergm_models/gof20.RDS")
rm(gof20)

## pooled prox ~ vi 95 day + vi 95 night + housedist*liveclose + gender match + age diff +gwesp(1, fixed) +gwnsp
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist*m$liveclose, attrname = "houseclose") + 
                   nodematch("gender") + absdiff("age" , pow=1) +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m21")
tmp <- data.frame(mod = "m21", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## pooled prox ~ vi 95 day + vi 95 night + housedist + gender match + age diff +gwesp(1, fixed) +gwnsp
mymod <- myergm((pooled_prox_new_binary_net ~ edges + 
                   edgecov(m$VI_95_day, attrname = "vi_95_day") + edgecov(m$VI_95_night, attrname = "vi_95_night") +
                   edgecov(m$housedist, attrname = "housedist") + 
                   nodematch("gender") + absdiff("age" , pow=1) +
                   gwesp(1, fixed = TRUE) + gwnsp()),
                "ergm_models/m22")
tmp <- data.frame(mod = "m22", aic = summary(mymod$model)$aic, 
                  accuracy = mymod$accuracy, specificity = mymod$specificity, sensitivity = mymod$sensitivity)
ergm_selection <- rbind(ergm_selection, tmp)

## check gof for model 
gof21 <- gof(readRDS("ergm_models/m21.RDS"))
gof21 
plot(gof21) 
saveRDS(gof21, "ergm_models/gof21.RDS")
rm(gof21)

## check gof for model 22
gof22 <- gof(readRDS("ergm_models/m22.RDS"))
gof22 
plot(gof22) 
saveRDS(gof22, "ergm_models/gof22.RDS")
rm(gof22)

## those look about the same so just go with model 21 (lower AIC)

## save the model selection table
ergm_selection <- ergm_selection %>% arrange(aic)
write.csv(ergm_selection, "ergm_selection_table.csv", row.names = F)

rm(mymod, tmp)

## check mcmc.diagnostic for best model with naming (14.1)
mcmc.diagnostics("ergm_models/m14.1.RDS")

## check mcmc.diagnostic for best model without naming (21)
mcmc.diagnostics(readRDS("ergm_models/m21.RDS"))

rm(ergm_selection, tmp, mymod)
########################################
## 2.8 simulated networks ####
########################################
m14.1 <- readRDS("ergm_models/m14.1.RDS")
#simulate 1000 networks from the model output. Note that the FULL observed prox network is the basis.
sim_m14.1 <-simulate(m14.1,
                     nsim = 1000, seed = 246,
                     burnin = 1000, interval = 1000,
                     basis = observed_prox_new_binary_net,
                     verbose = FALSE)
## save those
saveRDS(sim_m14.1, "ergm_models/simulated_networks_fullmodel.RDS")

## and simulate networks for model with no naming data
m21 <- readRDS("ergm_models/m21.RDS")
#simulate 1000 networks from the model output. Note that the FULL observed prox network is the basis.
sim_m21 <-simulate(m21,
                   nsim = 1000, seed = 246,
                   burnin = 1000, interval = 1000,
                   basis = observed_prox_new_binary_net,
                   verbose = FALSE)
## save those
saveRDS(sim_m21, "ergm_models/simulated_networks_nonamemodel.RDS")

rm(m14.1, m21)

########################################
## 2.9 derive and compare edge probabilities ####
########################################
# From these 1000 simulations, we can derive a "mean" network in which each edge is the mean of all 1000 
# binary values assigned to that dyad. Therefore, the mean network will be weighted, with edges ranging 
# from 0 to 1. Each weighted edge will represent the probability of close contact between a dyad. 

## for the best model (includes naming data)
simmatrix <- list()
for(i in 1:1000) {
  output <- as.matrix(sim_m14.1[[i]])
  simmatrix[[i]] <- output #create a list of 1000 matrices
}

# proportion of networks an edge exists in
edgeprobs <- Reduce("+", simmatrix)/1000
edgeprobs[lower.tri(edgeprobs)] <- NA
diag(edgeprobs) <- NA

## a couple of notable values
length(which(edgeprobs == 0)) #123
length(which(edgeprobs == (1/1000))) #28
length(which(edgeprobs == 1)) #487

range(edgeprobs, na.rm = T) #0 1
mean(edgeprobs, na.rm = T) #0.2739636
length(which(!is.na(edgeprobs)));choose(123,2)

## save those as an edgelist
edgeprobs_m14.1 <- reshape2::melt(edgeprobs) %>%
  rename(id1 = Var1, id2 = Var2, edgeprob = value) %>%
  drop_na()

write.csv(edgeprobs_m14.1, file = "ergm_mean_edge_prob_names.csv", row.names = F)

## for the best model (doesnt include naming data)
simmatrix <- list()
for(i in 1:1000) {
  output <- as.matrix(sim_m21[[i]])
  simmatrix[[i]] <- output #create a list of 1000 matrices
}

# proportion of networks an edge exists in
edgeprobs <- Reduce("+", simmatrix)/1000
edgeprobs[lower.tri(edgeprobs)] <- NA
diag(edgeprobs) <- NA

## a couple of notable values
length(which(edgeprobs == 0)) #125
length(which(edgeprobs == (1/1000))) #29
length(which(edgeprobs == 1)) #474

range(edgeprobs, na.rm = T) #0 1
mean(edgeprobs, na.rm = T) #0.2748245
length(which(!is.na(edgeprobs)));choose(123,2)

## save those as an edgelist
edgeprobs_m21 <- reshape2::melt(edgeprobs) %>%
  rename(id1 = Var1, id2 = Var2, edgeprob = value) %>%
  drop_na()

write.csv(edgeprobs_m21, file = "ergm_mean_edge_prob_nonames.csv", row.names = F)

## compare the 2
full_join(edgeprobs_m14.1, edgeprobs_m21, by = c("id1", "id2")) %>% 
  ggplot(aes(x = edgeprob.x, y = edgeprob.y)) +
  geom_point() + 
  geom_abline()
cor.test(edgeprobs_m14.1$edgeprob, edgeprobs_m21$edgeprob, method = "spearman") #rho = 0.97222 p < 2.2e-16

rm(edgeprob, output, simmatrix)

########################################
## 3. GLM ~ predict edge weights ####
########################################
require(tidyverse)
## pairs df
pairs <- read.csv("big_edgelist_Mandena1.csv")

## first check that the number of simultaneous points and prox are not correlated
ggplot(pairs, aes(x = n_simultaneous, y = prox)) +
  geom_point()+geom_smooth() +
  labs(title = "Check correlation n_simultaneous and prox", 
       subtitle = paste("Rho (spearmans) = ", round(cor(pairs$n_simultaneous, pairs$prox),4)))

completepairs <- pairs %>% 
  filter(n_simultaneous >240, !is.na(VI_95_cowear)) #this is half a day of gps points (20pt/hr*12hrs)
########################################
## 3.1 Data distributions and correlations   ####
########################################
## VI estimate spearman correlations  
tmp <- completepairs %>% filter(VI_95_cowear > 0)
plot(x = tmp$VI_95, tmp$VI_95_cowear)
plot(x = completepairs$VI_95, completepairs$VI_95_cowear)
plot(x = completepairs$prox, completepairs$VI_95_cowear)
# try spearman coef
cor.test(completepairs$VI_95, completepairs$VI_95_cowear, method = "spearman")
cor.test(completepairs$VI_50, completepairs$VI_50_cowear, method = "spearman")
cor.test(completepairs$prox, completepairs$VI_95, method = "spearman")
cor.test(completepairs$prox, completepairs$VI_95_cowear, method = "spearman")
cor.test(completepairs$prox, completepairs$VI_50, method = "spearman")
cor.test(completepairs$prox, completepairs$VI_50_cowear, method = "spearman")

## distribution of the variables  
completepairs %>% count(gender_match)
completepairs %>% count(liveclose) #1 means <= 25m apart
# completepairs %>% count(recip.name) 


ggplot(completepairs, aes(x = prox))+geom_histogram(binwidth = 0.001)+ggtitle("Proximity")
ggplot(completepairs[which(completepairs$prox != 0),], aes(x = prox))+geom_histogram(binwidth = 0.001)+ggtitle("Proximity > 0")
ggplot(completepairs, aes(x=kader:::cuberoot(VI_95), y=gender_match))+
  ggridges::geom_density_ridges()

tmp <- completepairs %>% 
  select(contains("95"), contains("50")) %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
ggplot(data = tmp, aes(x = variable, y = value))+geom_boxplot()+ggtitle("VI values")
ggplot(data = completepairs, aes(x = age_diff))+geom_boxplot()+ggtitle("Age difference")
# ggplot(data=completepairs[completepairs$n.named >0,], aes(x=as.factor(n.named)))+geom_bar()+ggtitle("n.named > 0")
ggplot(data=completepairs, aes(x=housedist))+geom_boxplot() +ggtitle("house distance")

########################################
## 3.2 train/test datasets   ####
########################################
set.seed(246)
train.index <- sample(1:nrow(completepairs), 0.8 * nrow(completepairs))
train.dat <- completepairs[train.index,]
test.dat <- completepairs[-train.index,]

########################################
## 3.3 model selection   ####
########################################

## 3.3.1 global model #####
## cowear/VI_95 and VI_50
co.xi <- tweedie::tweedie.profile(prox ~ kader:::cuberoot(VI_95_cowear) + kader:::cuberoot(VI_50_cowear)
                                  # +n.named +recip.name
                                  +gender_match +age_diff + housedist*liveclose, 
                                  data = completepairs, link.power = 0)
co <- glm(prox ~ kader:::cuberoot(VI_95_cowear) + kader:::cuberoot(VI_50_cowear) +n.named +recip.name+gender_match +age_diff + housedist*liveclose, data = completepairs, family = statmod::tweedie(var.power = co.xi$xi.max, link.power = 0))

## night day VI_95
nd.xi <- tweedie::tweedie.profile(prox ~ kader:::cuberoot(VI_95_day) + kader:::cuberoot(VI_50_night)
                                  # +n.named +recip.name
                                  +gender_match +age_diff + housedist*liveclose, 
                                  data = completepairs, link.power = 0)
nd <- glm(prox ~ kader:::cuberoot(VI_95_night) + kader:::cuberoot(VI_50_day) +n.named +recip.name+gender_match +age_diff + housedist*liveclose, data = completepairs, family = statmod::tweedie(var.power = nd.xi$xi.max, link.power = 0))

## 3.3.2 choose top model ####
### dredge 
options(na.action = "na.fail")
co.d <- MuMIn::dredge(co, rank = tweedie::AICtweedie)
nd.d <- MuMIn::dredge(nd, rank = tweedie::AICtweedie)
options(na.action = "na.omit")

### view the best models by AIC
subset(co.d, delta <2) 
subset(nd.d, delta <2) 

### choose model with delta <2 that minimizes degrees of freedom
### choose model with delta <2 that minimizes degrees of freedom
co_best <- MuMIn::get.models(co.d, 1)[[1]]
summary(co_best)
nd_best <- MuMIn::get.models(nd.d, 4)[[1]]
summary(nd_best)

saveRDS(list(model_co_best = co_best,
             model_nd_best = nd_best),
        "glm_best_models.RDS")

########################################
## 3.4 check residuals  ####
########################################
## vi cowear
qres <-statmod::qresid(co_best)
qqnorm(qres, las=1); qqline(qres) #as suggested by tweedie package author in "Generalized Linear Model with Examples in R" 2018, pg 469. ISBN 978-1-4419-0117-0
plot(qres ~ fitted(co_best), las=1, xlab="Fitted values", ylab="Quantile residuals" )
plot( cooks.distance(co_best), type="h", las=1, ylab="Cook's distance, D")
print("influential predictors")
q.inf <- influence.measures(co_best)
colSums(q.inf$is.inf)

## vi night/day
qres <-statmod::qresid(nd_best)
qqnorm(qres, las=1); qqline(qres) #as suggested by tweedie package author in "Generalized Linear Model with Examples in R" 2018, pg 469. ISBN 978-1-4419-0117-0
plot(qres ~ fitted(nd_best), las=1, xlab="Fitted values", ylab="Quantile residuals" )
plot( cooks.distance(nd_best), type="h", las=1, ylab="Cook's distance, D")
print("influential predictors")
q.inf <- influence.measures(nd_best)
colSums(q.inf$is.inf)

########################################
## 3.5 sens/spec/acc ####
########################################
## 3.5.1 function ####
ssa_func <- function(test.dat, train.dat, model, plot = TRUE){
  model <- update(model, data = train.dat)
  test.dat$predicted <- predict(model, test.dat, type = "response") # get the predictions
  
  r <- range(test.dat$predicted) #check the range is 0-1
  
  if (plot){
    plot(test.dat$predicted, test.dat$prox, xlab = "Predicted", ylab = "Observed");abline(0,1, col = "blue")
  }
  #if predictions >1 see which ones they are, then force them to become 1
  toobig <- test.dat[test.dat$predicted >1,]
  test.dat$predicted <- ifelse(test.dat$predicted>1, 1, test.dat$predicted)
  
  
  test.dat <- test.dat %>% mutate(plabels = ifelse(prox>0, 1,0)) #make a variable of the 'truth' variable
  pred <- ROCR::prediction(test.dat$predicted, test.dat$plabels)
  perf_sens <- ROCR::performance(pred, measure = "sens") #sensitivity specificity curves
  perf_spec <- ROCR::performance(pred, measure = "spec")
  if (plot) {
    plot(unlist(perf_sens@x.values), unlist(perf_sens@y.values), 
         type="l", lwd=2, ylab="Sensitivity", xlab="Cutoff", xlim = c(0,.02))
    par(new=TRUE)
    plot(unlist(perf_spec@x.values), unlist(perf_spec@y.values), 
         type="l", lwd=2, col='red', ylab="", xlab="", xlim = c(0,.02))
    axis(4, at=seq(0,1,0.2))
    mtext("Specificity",side=4, padj=-2, col='red')
    
  }
  sens = cbind(unlist(perf_sens@x.values), unlist(perf_sens@y.values))
  spec = cbind(unlist(perf_spec@x.values), unlist(perf_spec@y.values))
  th = sens[which.min(apply(sens, 1, function(x) min(colSums(abs(t(spec) - x))))), 1]
  
  test.dat <- test.dat %>% ## see where the threshold fall with real prox values
    mutate(threshold = as.factor(ifelse(prox <=th, 0, 1)))
  
  if (plot) {
    plot(test.dat$threshold, log(test.dat$predicted), xlab = "Prox < threshold", ylab = "log(predicted)")
    abline(h= log(th), col="blue")
  }
  
  zz = test.dat %>% filter(threshold == 0, predicted <=th) %>% nrow() #correct 0s
  zo = test.dat %>% filter(threshold == 0, predicted >th) %>% nrow() #incorrect 0s
  oz = test.dat %>% filter(threshold == 1, predicted <=th) %>% nrow() #incorrect 1s
  oo = test.dat %>% filter(threshold == 1, predicted >th) %>% nrow() #correct 1s
  ## ability to predict below threshold (0s) (with zi = ~1, zi =~VI)
  sens = oo/(oo+oz) #sensitivity 
  spec = zz/(zz+zo)# specificity
  ac = (zz+oo)/(zz+zo+oz+oo) #predictive accuracy
  
  ## rsq
  rss <- sum((test.dat$pred - test.dat$prox) ^2)
  tss <- sum((test.dat$prox - mean(test.dat$prox))^2)
  rsq <- 1 - rss/tss
  
  return(list(
    "model" = model$call,
    "range" = r,
    "threshold" = th,
    "sensitivity" = sens,
    "specificity" = spec, 
    "accuracy" = ac, 
    "Rsq" = rsq,
    "too_big_predictions" = toobig))
}

## 3.5.2 run function ####
ssa_co <- ssa_func(test.dat = test.dat, train.dat = train.dat, model = co_best, plot=TRUE)
sjPlot::tab_model(co_best) 

ssa_nd <- ssa_func(test.dat = test.dat, train.dat = train.dat, model = nd_best, plot=TRUE)
sjPlot::tab_model(nd_best) 

## 3.5.3 see if predicted values outside expected range (0-1) ####
ssa_co$too_big_predictions
ssa_nd$too_big_predictions

## 3.5.2 check sens/spec/acc on multiple iterations####
ssa <- data.frame()
for(i in 1:100){
  t.index <- sample(1:nrow(completepairs), 0.8 * nrow(completepairs))
  tr.dat <- completepairs[t.index,]
  te.dat <- completepairs[-t.index,]
  ssa_nd <- ssa_func(test.dat = te.dat, train.dat = tr.dat, model = nd_best, plot = FALSE)
  ssa_co <- ssa_func(test.dat = te.dat, train.dat = tr.dat, model = co_best, plot = FALSE)
  d <- data.frame(mod = c("nightday_GLM_tweedie", "cowear_GLM_tweedie"),
                  sens = c(ssa_nd$sensitivity, ssa_co$sensitivity), 
                  spec = c(ssa_nd$specificity, ssa_co$specificity), 
                  acc = c(ssa_nd$accuracy, ssa_co$accuracy), 
                  threshold = c(ssa_nd$threshold, ssa_co$threshold),
                  Rsq = c(ssa_nd$Rsq, ssa_co$Rsq))
  ssa <-rbind(ssa, d)
}

ssasum <- ssa %>% 
  group_by(mod) %>% 
  summarise(
    sensitivity = mean(sens), 
    specificity = mean(spec), 
    accuracy = mean(acc), 
    threshold = mean(threshold), 
    Rsq = mean(Rsq)) %>% 
  mutate(AIC = c(tweedie::AICtweedie(nd_best), tweedie::AICtweedie(co_best)))
ssasum ##choose cowear model

########################################
## 3.6 observed vs predicted plots  ####
########################################
completepairs$predicted <- predict(co_best, completepairs, type = "response")
completepairs$predicted <- ifelse(completepairs$predicted >1, 1, completepairs$predicted)
plot(completepairs$prox, completepairs$predicted, cex = 0.5,
     xlab = "Observed (log)", ylab = "Predicted (log)", log = "xy");abline(0,1, col = "blue")
title(main = "GLM tweedie Cowear",
      sub = paste( "Rho = ", round(cor(completepairs$prox, completepairs$predicted, 
                                       method = "spearman"),2), sep = " "))

########################################
## 3.8 expected proportion of 0s  ####
########################################
## This stems from the discussion and methods outlined in this thread:  
### <https://stats.stackexchange.com/questions/174121/can-a-model-for-non-negative-data-with-clumping-at-zeros-tweedie-glm-zero-infl>
test.dat <- completepairs
Phi <- summary(co_best)$dispersion  #dispersion parameter of model
Mu <- fitted(co_best)
Power <- co.xi$xi.max
Prob.Zero <- exp(-Mu^(2-Power) / Phi / (2-Power))
Prob.Zero <- sort(Prob.Zero)
mean(Prob.Zero) #mean proportion of predicted 0s
sd(Prob.Zero) #sd of proportion of predicted 0s
range(Prob.Zero) #range of proportion of predicted zeros
boxplot(Prob.Zero)

## figure out what the cutoff values would be there
test.dat$co_predicted <- predict(co_best, test.dat, type = "response")
pred <- sort(na.omit(test.dat$co_predicted))
cutoff <- round(mean(Prob.Zero)*length(pred))
## as compared to sens/spec function method
cf <-pred[cutoff] 
th <- ssasum[ssasum$mod == "cowear_GLM_tweedie",]$threshold
round(length(pred[pred<th])/length(pred),4)

########################################
## 3.9 predict close contact edge weights ####
########################################
## use cowear models used to extrapolate to the whole VIs to predict edge weights
#change the column names so model will work
pred <- pairs %>% select(-VI_95_cowear, -VI_50_cowear) %>% rename(VI_95_cowear = VI_95, VI_50_cowear = VI_50) 
## GLM cowear
pred$edgeweight <- predict(co_best, pred, type = "response")
## if predict value > 1 change to 1
pred$edgeweight <- ifelse(pred$edgeweight >1, 1, pred$edgeweight)

##save it
pred <- pred %>% 
  select(id1, id2, edgeweight)
write.csv(pred, file = "glm_edgeweights.csv", row.names = F)

########################################
## 4. Close contact network ####
########################################
woregps <- read.csv(file = "woregps.csv")[,1]
## mean edge probabilities
edgeprobs <- read.csv(file = "ergm_mean_edge_prob_nonames.csv")
## predicted edge weights
pred <- read.csv(file = "glm_edgeweights_noname.csv")
## make those into 1 df
cc.edges <- full_join(edgeprobs, pred, by = c("id1", "id2")) %>% 
  filter(id1 %in% woregps, id2 %in% woregps) 
rm(edgeprobs, pred)

## simulated networks
sim <- readRDS(file = "ergm_models/simulated_networks_nonamemodel.RDS")

########################################
## 4.1 compare edges ####
########################################
## Compare the ERGM imputed edge probabilities to GLM imputed prox values
ggplot(cc.edges, aes(x = edgeweight, y = edgeprob)) +
  geom_point() +
  xlab("GLM Edge Weight") + ylab("ERGM Edge Probability") + 
  xlim(0,0.05) +
  geom_vline(xintercept = 0.002 , linetype = "dashed")

## spearman correlations
round(cor(cc.edges[,c("edgeweight", "edgeprob")], method = "spearman", use = "complete.obs"), 2)

########################################
## 4.1 combine ergm and glm outputs ####
########################################
# We can combine both the ERGM and GLM imputed vales to obtain a more robust simulated network of close contact. 
# The ERGM outputs "probability" of close contact (from 0 to 1) while the GLM outputs contact rates 
# (prox, also from 0 to 1). In addition, we obtain 1000 binary networks from the ERGM, and one weighted network
# from the GLM. Therefore, to combine them, we multiply each binary ERGM network by the GLM network and repeat 
# this 1000 times. Then, some edges that have an ERGM probability of 1 will be assigned their maximum GLM edge weight, 
# while edges with probability < 1 will be downweighted. 

## the edge weight * the probability an edge exists ####
cc.edges <- cc.edges %>% 
  mutate(imputed.weight = edgeweight * edgeprob)
write.csv(cc.edges, file = "close_contact_edges_nonames.csv", row.names = F)

## the edge weight on each of the simulated networks  ####
## first make the glm predicted weights into a matrix
cc.tmp <- cc.edges %>% 
  rename(id1 = id2, id2 = id1) %>% 
  rbind(cc.edges)
weight.matrix <- reshape2::dcast(cc.tmp, id1~id2, value.var = "edgeweight") 
row.names(weight.matrix) <- weight.matrix$id1
weight.matrix <- weight.matrix %>% 
  dplyr::select(-id1) %>% 
  as.matrix
rm(cc.tmp)

## make upper and lower symetric to match simulated networks format
which(upper.tri(weight.matrix) != t(lower.tri(weight.matrix))) #to check same
# weight.matrix[lower.tri(weight.matrix)] <- t(weight.matrix)[lower.tri(weight.matrix)]

## make the network.list into list of matrices again and multiple each by predicted weight
cc.matrix <- list()
for(i in 1:length(sim)) {
  output <- as.matrix(sim[[i]])
  cc.matrix[[i]] <- output * weight.matrix
}

########################################
## 4.2 make graphs ####
########################################
## calculate metric on each of the simulated matrices to get a distribution of values

## make a list of graphs
igraph_list <- list()

for(i in 1:length(sim)){
  ## make a weighted graph object
  g <- igraph::graph_from_adjacency_matrix(cc.matrix[[i]], mode = "undirected", 
                                           weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = TRUE) #ICHAGNED DIAG TO f TO MATCH CHANGE ABOVE AND DON'T NEED SIMIPLIFY
  ## simplify to be safe
  g <- igraph::simplify(g)
  ## add to list
  igraph_list[[i]] <- g
}

########################################
## 4.3 network wide stats ####
########################################
widestats <- data.frame()

for (i in 1:length(sim)) {
  df <- data.frame(
    connected_strong = igraph::is.connected(igraph_list[[i]], mode = "strong"),
    connected_weak = igraph::is.connected(igraph_list[[i]], mode = "weak"),
    diameter = igraph::diameter(igraph_list[[i]]),
    avg_distance = igraph::average.path.length(igraph_list[[i]]),
    density = igraph::graph.density(igraph_list[[i]]),
    transitivity_global = igraph::transitivity(igraph_list[[i]], type = "global"),
    louvain_modularity = igraph::modularity(igraph::cluster_louvain(igraph_list[[i]])))
  
  widestats <- rbind(widestats, df)
}

widestats <- data.frame(
  connected_strong = mean(widestats$connected_strong),
  connected_weak = mean(widestats$connected_weak),
  diameter = mean(widestats$diameter),
  avg_distance = mean(widestats$avg_distance),
  density = mean(widestats$density),
  transitivity_global = mean(widestats$transitivity_global),
  louvain_modularity = mean(widestats$louvain_modularity),
  
  sd_diameter = sd(widestats$diameter),
  sd_avg_distance = sd(widestats$avg_distance),
  sd_density = sd(widestats$density),
  sd_transitivity_global = sd(widestats$transitivity_global),
  sd_louvain_modularity = sd(widestats$louvain_modularity)
)

## save it
write.csv(widestats, file = "widestats_close_contact_nonames.csv", row.names = F)
rm(widestats)

## number of components in each graph
comp <- c()

for (i in 1:length(sim)) {
  clu <- igraph::components(igraph_list[[i]])
  n <- length(igraph::groups(clu))
  comp <- c(comp, n)
}

table(comp) 
rm(comp, clu, n)

########################################
## 4.4 Strength centrality ####
########################################

node.strength <- data.frame()

for (i in 1:length(sim)) {
  s <- igraph::strength(igraph_list[[i]])
  s <- as.data.frame(s) %>% 
    mutate(id = rownames(.))
  node.strength <- rbind(node.strength, s)
}

## to plot someone
node.strength %>% 
  filter(id == unique(node.strength$id)[1]) %>% #change the 1 to plot someone else
  ggplot(aes(x = s)) +
  geom_density()

## boxplots of everyone
ggplot(node.strength, aes(x = reorder(id, s, FUN = median), y = s)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme(axis.text.y=element_blank())

## coef of variation for each individual
node.strength.var <- node.strength %>% 
  group_by(id) %>% 
  summarise(mean = mean(s), sd = sd(s), 
            median = median(s), 
            norm = shapiro.test(s)$p.value) %>% 
  mutate(cv = sd/mean)

ggplot(node.strength.var, aes(x = cv)) + geom_histogram(binwidth = .01 , color = "gray" , fill = "blue")
ggplot(node.strength.var, aes(x = median)) + geom_histogram(binwidth = .01 , color = "gray" , fill = "blue")

########################################
## 4.5 Degree centrality ####
########################################
node.degree <- data.frame()

for (i in 1:length(sim)) {
  d <- igraph::degree(igraph_list[[i]])
  d <- as.data.frame(d) %>% 
    mutate(id = rownames(.))
  node.degree <- rbind(node.degree, d)
}

## to plot someone 
node.degree %>% 
  filter(id == unique(node.degree$id)[1]) %>% #change the 1 to plot someone else
  ggplot(aes(x = d)) +
  geom_density()

## boxplots of everyone
ggplot(node.degree, aes(x = reorder(id, d, FUN = median), y = d)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme(axis.text.y=element_blank())

########################################
## 4.6 Eigenvector centrality ####
########################################
node.eigen <- data.frame()

for (i in 1:length(sim)) {
  # eigenvector not well defined for disconnected graphs so calculate for largest component
  sub_graph <- igraph::decompose(igraph_list[[i]], mode = "weak", max.comps = 1, min.vertices = 5)[[1]]
  e <- igraph::eigen_centrality(sub_graph)$vector
  e <- as.data.frame(e) %>% 
    mutate(id = rownames(.))
  node.eigen <- rbind(node.eigen, e)
}

## to plot someone 
node.eigen %>% 
  filter(id == unique(node.eigen$id)[1]) %>% #change the 1 to plot someone else
  ggplot(aes(x = e)) +
  geom_density()

## boxplots of everyone
ggplot(node.eigen, aes(x = reorder(id, e, FUN = median), y = e)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme(axis.text.y=element_blank())

########################################
## 4.7 Betweenness centrality ####
########################################
node.btwn <- data.frame()

for (i in 1:length(sim)) {
  # betweenness is well defined for disconnected graphs so don't need to calculate for largest component
  # COURTNEY WHY?
  sub_graph <- igraph::decompose(igraph_list[[i]], mode = "weak", max.comps = 1, min.vertices = 5)[[1]]
  b <- igraph::betweenness(sub_graph, weights = (1 - (igraph::E(sub_graph)$weight)/1.01)) #why divide by 1.01?
  ## node betweeeness in igraph views weights as distances so want the inverse
  b <- as.data.frame(b) %>% 
    mutate(id = rownames(.))
  node.btwn <- rbind(node.btwn, b)
}

## to plot someone 
node.btwn %>% 
  filter(id == unique(node.btwn$id)[1]) %>% #change the 1 to plot someone else
  ggplot(aes(x = b)) +
  geom_density()

## boxplots of everyone
ggplot(node.btwn, aes(x = reorder(id, b, FUN = median), y = b)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme(axis.text.y=element_blank())

########################################
## 4.8 Median centrality df ####
########################################
## make dfs of medians for each node
node.strength <- node.strength %>% 
  group_by(id) %>% 
  summarise(med_strength = median(s))
node.degree <- node.degree %>% 
  group_by(id) %>% 
  summarise(med_deg = median(d))
node.eigen <- node.eigen %>% 
  group_by(id) %>% 
  summarise(med_eig = median(e))
node.btwn <- node.btwn %>% 
  group_by(id) %>% 
  summarise(med_bet = median(b))

## join together
node.centrality <- node.strength %>% 
  left_join(node.degree, by = "id") %>% 
  left_join(node.eigen, by = "id") %>% 
  left_join(node.btwn, by = "id")

## save it
write.csv(node.centrality, file = "centrality_closecontact_nonames.csv", row.names = F)
