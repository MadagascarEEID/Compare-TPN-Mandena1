## Clean tracking data

require(tidyverse)
require(sp)
require(sf)
require(lubridate)
require(SDLfilter)


## first function flags fridays, duplicates, loops and
## points based on speed, timelag between fixes
## second function flags points outside study area and days not wearing GPS

gps_flags_1 <- function(raw, timelag = 0, prj, tmax = 30, smax = 12, altmin = 0, altmax = 2000) { 
  #timelag = 0 or 3 min to be considered duplicate
  #tmax = maximum amount of time in min t-1 to t to allow into dataset (less accurate when large gap) (30)
  #smax = maximum speed km/hr lead and lag average to allow a point to have (12)
  # altmax/altmin are the min/max alt to allow in 
  ###### Camp Marojejia at 775 m also looks reasonable on G Earth
  col1 <-colnames(raw)[1]
  
  df <- raw %>% rename(id = colnames(raw)[1]) %>% 
    mutate(id_day = paste(id, yday(raw$timestamp), sep="_")) %>% 
    mutate(id_week = paste(id, gps_week, sep = "_")) %>% 
    ## flag Fridays (GPS distribution day, don't want to include)
    ## devices only recorded on Fridays the first few weeks of study then settings were changed
    mutate(flag_fri = if_else(wday(raw$timestamp, label = TRUE) == "Fri", TRUE, FALSE))%>% 
    mutate(qi = row_number()+10000000) #unique row id
  
  ## SDL filter for duplicates (dup) and loops (dd)
  sdl <- df %>% 
    arrange(id_week, timestamp) %>% 
    dplyr::select(id_week, timestamp, latitude, longitude, qi) %>% 
    rename(id = id_week, DateTime = timestamp, lat = latitude, lon = longitude)
  
  
  ## loop through every id (person-day) ##or week...
  names <- unique(sdl$id)
  clean <- data.frame()
  
  ## filter
  for (i in 1:length(names)) { ##changing this to parallel would help speed a lot
    
    tmp <- filter(sdl, id %in% names[i]) #subset sdl df
    dup <- dupfilter(tmp, step.time = timelag/60, step.dist = 0, conditional = TRUE) #could use dplyr distinct instead...
    V <- vmax(dup, prob = .99) # set max one-way linear speed as upper 1%
    VLP <- suppressWarnings(vmaxlp(dup, prob = .95)) # set max loop speed as upper 5%
    dd <- suppressWarnings(ddfilter(dup, vmax=V, vmaxlp=VLP)) # filter points based on ests
    
    ## note what was filtered so don't throw it away
    dup$flag_dup <- FALSE
    dd$flag_dd <- FALSE
    wf <- tmp %>% 
      left_join(dd) %>% 
      left_join(dup) %>% 
      mutate(flag_dd = if_else(is.na(flag_dd), TRUE, FALSE)) %>% 
      mutate(flag_dup = if_else(is.na(flag_dup), TRUE, FALSE)) 
    
    ## save to new df
    clean <- rbind(clean, wf)
  }
  
  gps <- clean %>%
    dplyr::select(qi, flag_dd, flag_dup) %>% 
    mutate(qi = as.integer(qi)) %>% 
    # rename() %>% 
    dplyr::full_join(df) %>% 
    dplyr::select(-qi)
  
  ## filter again based on own settings
  
  ## make spatial
  gps <- st_as_sf(gps, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84" )
  gps <- st_transform(gps, prj) 
  ## Get the coordinates into a column
  gps[,c("xcoord","ycoord")] <- st_coordinates(gps)
  
  
  ## flag points based on own filters
  gpsClean <- gps %>% 
    arrange(id, timestamp) %>%  # Make sure everything is organized
    group_by(id_week) %>% # work within the same group
    mutate(dist1 = sqrt((lead(xcoord,default = NA) - xcoord)^2 + (lead(ycoord,default = NA) - ycoord)^2)) %>% # I don't trust gDistance, this is just to get distance between consecutive points manually
    mutate(dist2 = sqrt((lag(xcoord,default = NA) - xcoord)^2 + (lag(ycoord,default = NA) - ycoord)^2)) %>%
    mutate(difft1 = as.numeric(difftime(lead(timestamp,default = NA),timestamp,units = "hours"))) %>%
    mutate(difft2 = as.numeric(difftime(timestamp,lag(timestamp,default = NA),units = "hours"))) %>%
    mutate(speed = (((dist1)/difft1)+((dist2)/difft2))/2) %>% ## speed in m/hour
    mutate(flag_badpt = if_else(difft1 > tmax/60 | #if >tmax min between points
                                  # speed > smax*1000 | # speed cutoff should have been done by ddfilter
                                  altitude> altmax | 
                                  altitude< altmin | # Mandena = 82 m lowest local point 60 m
                                  # consider adding angle filter
                                  flag_dd == TRUE |
                                  flag_dup == TRUE,
                                TRUE, FALSE, missing = FALSE)) %>%
    mutate(flag_still = if_else(speed < 500, TRUE, FALSE, missing = FALSE)) %>% #if moving less than .5km/hr
    as.data.frame() %>% 
    dplyr::select(id, id_week, id_day, gps_week, timestamp, xcoord, ycoord, altitude, dist1, dist2, difft1, difft2, speed, 
                  flag_fri, flag_still, flag_badpt, flag_dd, flag_dup, geometry)
  colnames(gpsClean)[1] <- col1
  
  
  gpsClean %>% 
    dplyr::select(contains("flag")) %>% 
    summary() %>% 
    print()
  return(gpsClean)
}


require(adehabitatHR)
## second filter based on grid and wearing GPS #####

gps_flags_2 <- function(flagsoutput, type, grid, prj){
  
  tmp <- flagsoutput #%>% dplyr::select(-geometry)
  
  df <- tmp %>% 
    dplyr::filter(flag_badpt == FALSE,
                  flag_fri==FALSE)
  
  ## make spatial
  coordinates(df) <- ~xcoord + ycoord
  raster::crs(df) <- prj
  
  ## select points within the grid area
  gps <- as.data.frame(df[grid,])
  gps$flag_grid <- FALSE
  
  
  df <- as.data.frame(df) %>% 
    full_join(gps) %>% 
    mutate(flag_grid = if_else(is.na(flag_grid), TRUE, FALSE))
  if (type=="animal") {
    df.final <- df %>% 
      full_join(tmp) %>% 
      dplyr::select(-flag_dup,-flag_dd) #both captured in flag_badpt
    return(df.final)
  }
  if(type == "person"){
    ## flag days where the person traversed <1ha
    day.area <- df %>% 
      dplyr::group_by(id_day) %>% dplyr::filter(n()>5)
    coordinates(day.area) <- ~xcoord + ycoord
    raster::crs(day.area) <- prj
    day.filter <- mcp.area(day.area[,"id_day"], percent = 99, plotit = F) #could change or look at multiple %
    day.filter <- day.filter %>% 
      rownames_to_column("perc") %>%
      pivot_longer(cols = -c("perc"), names_to = "id", values_to = "size") %>%
      dplyr::filter(perc == "99" & size >1.) %>%  # if change above also change here
      pull(id)
    df.day <- df %>% 
      mutate(flag_day = if_else(id_day %in% day.filter, FALSE, TRUE))
    
    df.final <- df.day %>% 
      # dplyr::select(-geometry) %>% 
      full_join(tmp) %>% 
      dplyr::select(-flag_dup,-flag_dd) #both captured in flag_badpt
    return(df.final)
  }
}

# ## plotting function
# myplot <- function(df, id, which.filter, prj){ #filter=1 or 2, id = 1, 2, 3...
#   df <- df %>% rename(who = names(df)[1]) #so animal and people have same name
# 
#   x <- df[df$who == unique(df$who)[id],]
# 
#   if (which.filter=="gps_flags_1") {
#     spdf <- SpatialPointsDataFrame(x[,c("xcoord", "ycoord")],
#                                    data = x[,!names(x) %in% c("xcoord", "ycoord")],
#                                    proj4string = CRS(prj))
#     plot(spdf, pch=20, cex=.75, col="black")
#     plot(spdf$geometry[spdf$flag_fri==TRUE], pch=20, cex=.5, col="grey", add=T)
#     plot(spdf$geometry[spdf$flag_badpt==TRUE], pch=20, cex=.5, col="red", add=T)
#     plot(spdf$geometry[spdf$flag_dup==TRUE], pch=20, cex=.5, col="orange", add=T)
#     plot(spdf$geometry[spdf$flag_dd==TRUE], pch=20, cex=.5, col="yellow", add=T)
#     legend("topright", legend = c("Friday", "badpt", "dup", "dd"),
#            col = c("grey", "red", "orange", "yellow"), pch=20)
#     title(unique(df$who)[id])
#   }
#   if (which.filter=="gps_flags_2") {
#     # x <- x %>% dplyr::filter(flag_fri==FALSE, flag_badpt==FALSE)
#     spdf <- SpatialPointsDataFrame(x[,c("xcoord", "ycoord")],
#                                    data = x[,!names(x) %in% c("xcoord", "ycoord")],
#                                    proj4string = CRS(prj))
#     plot(spdf, pch=20, cex=.75, col="black")
#     plot(spdf$geometry[spdf$flag_fri==TRUE], pch=20, cex=.5, col="grey", add=T)
#     plot(spdf$geometry[spdf$flag_badpt==TRUE], pch=20, cex=.5, col="red", add=T)
#     plot(spdf$geometry[spdf$flag_grid==TRUE], pch=20, cex=.5, col="green", add=T)
#     plot(spdf$geometry[spdf$flag_day==TRUE], pch=20, cex=.5, col="purple", add=T)
#     legend("topright", legend = c("Friday", "badpt", "outside grid", "day not worn"),
#            col = c("grey", "red", "green", "purple"), pch=20)
#     title(unique(df$who)[id])
#   }
# }
