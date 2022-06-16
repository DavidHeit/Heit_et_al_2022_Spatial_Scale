 
#=========================================================================================#
#  THE SPATIAL SCALING AND INDIVIDUALITY OF HABITAT SELECTION IN A WIDESPREAD UNGULATE    #
#=========================================================================================#
#                               ECOGRAPHY, DATE, VOL, ISS                                 #

                                     # CODE AUTHOR #
                      # DAVID R HEIT - UNIVERSITY OF NEW HAMPSHIRE #

                                    # PAPER AUTHORS #
 # DAVID R HEIT, JOSHUA J. MILLSPAUGH, JON T. MCROBERTS, KEVYN H. WISKIRCHEN, JASON A. SUMNERS, 
# JASON L. ISABELLE, BARBARA J. KELLER, AARON M. HILDRETH, ROBERT A. MONTGOMERY, REMINGTON J. MOLL

 

#### Libraries          ##  CITATIONS
library(raster)         ##  Hijmans et al. 2013
library(tidyverse)      ##  Wickham 2017
library(lubridate)      ##  Grolemun & Wickham 2011
library(amt)            ##  Signer et al. 2019
library(sf)             ##  Pebesma 2018
library(MuMIn)          ##  Barton 2020
library(qpcR)           ##  Spiess 2018
library(rgdal)          ##  Bivand et al. 2020
library(ctmm)           ##  Fleming & Calabrese 2021
library(adehabitatHR)   ##  Calenge 2006
library(maptools)       ##  Bivand & Lewin-Koh 2021
library(spatstat)       ##  Baddeley et al. 2014
library(abind)          ##  Plate & Heiberger 2016
library(viridis)        ##  Garnier 2018
library(usdm)           ##  Naimi et al. 2014
library(zoo)            ##  Zeileis & Grothendieck 2005
library(ggpubr)         ##  Kassambara 2020
library(adehabitatLT)   ##  Calenge 2006
library(rgeos)          ##  Bivand & Rundel 2021

#======================================================================#
####   C R E A T E   S C A L E D    H A B I T A T   R A S T E R S   ####
#======================================================================#
## Make list of scale factors for the moving windows. These numbers
# correspond to the side lengths (in cells) for the moving window. So 3 would
# indicate a 3x3 window, and 166 would indicate a 166 x 166 window. 
cells <- seq(1,129,by=2)   # range of scales to create (65 total)
res   <- 30                # resolution of base raster (m)

#---------------------------#
#  F O R E S T   C O V E R  #
#---------------------------#
forest.n <- raster("Forest_North.tif") # 30m raster for forest cover in N study site
forest.s <- raster("Forest_South.tif") # 30m raster for forest cover in S and SE sites
 

for (i in 1:length(cells)){
  print(cells[i]*res)
  print(Sys.time())
  
  # create the scaled raster
  r <- focal(forest.n, w = matrix(1/cells[i]^2, 
                                  nc = cells[i], 
                                  nr = cells[i])) 
  writeRaster(r,paste("C:/Users/dheit/OneDrive/Documents/Missouri_Rasters/Forest_North_", ## output directory
                      cells[i]*res,".tif",sep = ""), overwrite=T)
  
  r <- focal(forest.s, w = matrix(1/cells[i]^2, 
                                  nc = cells[i], 
                                  nr = cells[i])) 
  writeRaster(r,paste("C:/Users/dheit/OneDrive/Documents/Missouri_Rasters/Forest_South_", ## output directory
                      cells[i]*res,".tif",sep = ""), overwrite=T)
} 

#--------------------------------#
#    R O A D   D I S T A N C E   #
#--------------------------------#
## Load starting rasters
roads.n.dist <- raster("roads_n_dist.tif") # 30m raster for road distance in N study site
roads.s.dist <- raster("roads_s_dist.tif") # 30m raster for road distance in S and SE sites

### Road Distance
for (i in 1:length(cells)){
  
  # create the scaled raster
  print(paste0("North All ",cells[i]*res," ", Sys.time()))
  r <- focal(roads.n.dist, w = matrix(1/cells[i]^2, nc = cells[i], nr = cells[i])) 
  writeRaster(r,paste("C:/Users/dheit/OneDrive/Documents/Missouri_Rasters/Roads_North_", ## output directory
                      cells[i]*30,".tif",sep = ""), overwrite=T)
  
  print(paste0("South All ",cells[i]*res, " ", Sys.time(), sep = " "))
  r <- focal(roads.s.dist, w = matrix(1/cells[i]^2, nc = cells[i], nr = cells[i])) 
  writeRaster(r,paste("C:/Users/dheit/OneDrive/Documents/Missouri_Rasters/Roads_South_", ## output directory
                      cells[i]*30,".tif",sep = ""), overwrite=T)
} 




#====================================#
####      L O A D   D A T A       ####
#====================================#

### Read-in Data
deer.all <- read.csv("Heit_2022_ECOG_deer_data.csv")
## GPS locations of deer locations
## Required columns: 
# id (individual), x (longitude), y (latitude), 
#  t (date/time), other covariates can be added after these

## Change Ids to factors
deer.all$id <- factor(deer.all$id)

#========================================#
####        P R E P  D A T A          ####
#========================================#

## Coerce dates to proper format
deer.all$t <- ymd_hms(deer.all$t)

## To increase efficiency we cropped habitat rasters as small as possible
#  This left us with two raster extents - one for the north and one for both the
#  south and southeast sites.

deer.n <- deer.all[deer.all$site=="North",]
deer.s <- deer.all[deer.all$site=="South",] # Southeast site included with South

## Make track and resample to 5 hour intervals and simulate
# random steps (n = 5)

# North
deer.trk.n <- deer.n %>% drop_na() %>%
  make_track(.x=x,.y=y,.t=t,id=id_yr_sn, crs = CRS("+init=epsg:5070"), all_cols = T) %>% nest(-id_yr_sn) %>%
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = hours(5), tolerance = minutes(15)) %>% steps_by_burst() %>%
      random_steps(n = 5)))

# South
deer.trk.s <- deer.s %>% drop_na() %>%
  make_track(.x=x,.y=y,.t=t,id=id_yr_sn, crs = CRS("+init=epsg:5070"), all_cols = T) %>% nest(-id_yr_sn) %>%
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = hours(5), tolerance = minutes(15)) %>% steps_by_burst() %>%
      random_steps(n = 5))) 


#============================================#
####          P H A S E  O N E            ####
#--------------------------------------------#
###     I N T E G R A T E D   S S F        ###
#============================================#

# prep data for loop
cells <- seq(1,129,by=2) ## 65 scales

#--------------------------------------#
###     F O R E S T   C O V E R      ###
#--------------------------------------#
### NORTH - Forest --------------------------------####
# results holder array
results_for_n <- array(NA, dim = c(nrow(deer.trk.n), # individuals
                                   length(cells), # scales
                                   9)) # parameters
for (i in 1:length(cells)){
  print(cells[i]*30)
  
  # create the scaled raster
  forest <- raster(paste("C:\\Users\\dh1157\\Documents\\Missouri_Rasters\\Forest_North_",cells[i]*30,".tif",sep = ""))
  names(forest) <- "forest"
  
  # extract covariates and fit issf
  mod <- deer.trk.n %>% 
    mutate(ssf = lapply(steps, function(x){
      x %>% 
        extract_covariates(forest) %>% 
        mutate(sl_2 = ifelse(sl_ > 0, sl_, 0.001)) %>% 
        fit_issf(case_ ~ scale(forest)[,1] + log(abs(sl_2)) + strata(step_id_)) # 'road' name created above
    }))
  
  # store results: rows = individuals, columns = parameter values, 3rd dim = parameter types
  results_for_n[,i,1] <- unlist(lapply(mod$ssf,AIC))                                # AIC
  results_for_n[,i,2] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,1]})) # coef  ; row 1 is for forest
  results_for_n[,i,3] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,2]})) # exp(coef)
  results_for_n[,i,4] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,3]})) # se(coef)
  results_for_n[,i,5] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,5]})) # p-val
  results_for_n[,i,6] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,1]})) # coef  ; row 2 is for log sl
  results_for_n[,i,7] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,2]})) # exp(coef)
  results_for_n[,i,8] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,3]})) # se(coef)
  results_for_n[,i,9] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,5]})) # p-val
  rm(forest)
  rm(mod)
  print(Sys.time())
}


#### South - Forest -----------------------------------####
results_for_s <- array(NA, dim = c(nrow(deer.trk.s), # individuals
                                   length(cells), # scales
                                   9)) # parameters
for (i in 1:length(cells)){
  print(cells[i]*30)
  
  # create the scaled raster
  forest <- raster(paste("C:\\Users\\dh1157\\Documents\\Missouri_Rasters\\Forest_South_",cells[i]*30,".tif",sep = ""))
  names(forest) <- "forest"
  
  # extract covariates and fit issf
  mod <- deer.trk.s %>% 
    mutate(ssf = lapply(steps, function(x){
      x %>% 
        extract_covariates(forest) %>% 
        mutate(sl_2 = ifelse(sl_ > 0, sl_, 0.001)) %>% 
        fit_issf(case_ ~ scale(forest)[,1] + log(abs(sl_2)) + strata(step_id_)) # 'road' name created above
    }))
  
  # store results: rows = individuals, columns = parameter values, 3rd dim = parameter types
  results_for_s[,i,1] <- unlist(lapply(mod$ssf,AIC))                                # AIC
  results_for_s[,i,2] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,1]})) # coef  ; row 1 is for forest
  results_for_s[,i,3] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,2]})) # exp(coef)
  results_for_s[,i,4] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,3]})) # se(coef)
  results_for_s[,i,5] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,5]})) # p-val
  results_for_s[,i,6] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,1]})) # coef  ; row 2 is for log sl
  results_for_s[,i,7] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,2]})) # exp(coef)
  results_for_s[,i,8] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,3]})) # se(coef)
  results_for_s[,i,9] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,5]})) # p-val
  rm(forest)
  rm(mod)
  print(Sys.time())
}


#----------------------------------------------#
#### North - All roads --------------------------------####
results_road_all_n <- array(NA, dim = c(nrow(deer.trk.n), # individuals
                                   length(cells), # scales
                                   9)) # parameters
for (i in 1:length(cells)){
  print(cells[i]*30)
  
  # create the scaled raster
  road <- raster(paste("C:\\Users\\dh1157\\Documents\\Missouri_Rasters\\Roads_North_",cells[i]*30,".tif",sep = ""))
  names(road) <- "road_all"
  
  # extract covariates and fit issf
  mod <- deer.trk.n %>% 
    mutate(ssf = lapply(steps, function(x){
      x %>% 
        extract_covariates(road) %>% 
        mutate(sl_2 = ifelse(sl_ > 0, sl_, 0.001)) %>% 
        fit_issf(case_ ~ log(abs(scale(road_all)[,1])) + log(abs(sl_2)) + strata(step_id_)) # 'road' name created above
    }))
  
  # store results: rows = individuals, columns = parameter values, 3rd dim = parameter types
  results_road_all_n[,i,1] <- unlist(lapply(mod$ssf,AIC))                                # AIC
  results_road_all_n[,i,2] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,1]})) # coef  ; row 1 is for forest
  results_road_all_n[,i,3] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,2]})) # exp(coef)
  results_road_all_n[,i,4] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,3]})) # se(coef)
  results_road_all_n[,i,5] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,5]})) # p-val
  results_road_all_n[,i,6] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,1]})) # coef  ; row 2 is for log sl
  results_road_all_n[,i,7] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,2]})) # exp(coef)
  results_road_all_n[,i,8] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,3]})) # se(coef)
  results_road_all_n[,i,9] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,5]})) # p-val
  rm(road)
  rm(mod)
  print(Sys.time())
}

#### South - All roads --------------------------------------####
results_road_all_s <- array(NA, dim = c(nrow(deer.trk.s), # individuals
                                        length(cells), # scales
                                        9)) # parameters
for (i in 1:length(cells)){
  print(cells[i]*30)
  
  # create the scaled raster
  road <- raster(paste("C:\\Users\\dh1157\\Documents\\Missouri_Rasters\\Roads_South_",cells[i]*30,".tif",sep = ""))
  names(road) <- "roads"
  
  # extract covariates and fit issf
  mod <- deer.trk.s %>% 
    mutate(ssf = lapply(steps, function(x){
      x %>% 
        extract_covariates(road) %>% 
        mutate(sl_2 = ifelse(sl_ > 0, sl_, 0.001)) %>% 
        mutate(roads = ifelse(roads > 0, roads, 0.001)) %>% 
        fit_issf(case_ ~ scale(log(roads))[,1] + log(sl_2) + strata(step_id_)) # 'road' name created above
    }))
  
  # store results: rows = individuals, columns = parameter values, 3rd dim = parameter types
  results_road_all_s[,i,1] <- unlist(lapply(mod$ssf,AIC))                                # AIC
  results_road_all_s[,i,2] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,1]})) # coef  ; row 1 is for forest
  results_road_all_s[,i,3] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,2]})) # exp(coef)
  results_road_all_s[,i,4] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,3]})) # se(coef)
  results_road_all_s[,i,5] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[1,5]})) # p-val
  results_road_all_s[,i,6] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,1]})) # coef  ; row 2 is for log sl
  results_road_all_s[,i,7] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,2]})) # exp(coef)
  results_road_all_s[,i,8] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,3]})) # se(coef)
  results_road_all_s[,i,9] <- unlist(lapply(mod$ssf,function(x){coef(summary(x))[2,5]})) # p-val
  rm(road)
  rm(mod)
  print(Sys.time())
}

#### Save results  ---------------------------------####
saveRDS(results_for_n,"results_forest_north.rds")
saveRDS(results_for_s,"results_forest_south.rds")
saveRDS(results_road_all_n,"results_road_all_north.rds")
saveRDS(results_road_all_s,"results_road_all_south.rds")




#============================================#
####        P H A S E    T W O            ####
#--------------------------------------------#
###      L I N E A R    M O D E L S        ###
#============================================#

#### Read in data
results_for_n <- readRDS("results_forest_north.rds")
results_for_s <- readRDS("results_forest_south.rds")
results_road_all_n <- readRDS("results_road_all_north.rds")
results_road_all_s <- readRDS("results_road_all_south.rds")


# REFERENCE:  results[ individual , scale , parameter]

## Data frame of Scales chosen by AIC ----------------------------#
## with individual-level covariates added back in

scale.road.df <- data.frame()
cells <- seq(1,129,by=2)
scales <- cells*30

## Roads North--------------#
deer.names <- unique(deer.n$id_yr_sn)
for (i in 1:length(deer.names)) {
  num.r <- which.min(results_road_all_n[i,,1])
  deer.i <- deer.n[deer.n$id_yr_sn==deer.names[i],]
  df <- data.frame(id = deer.names[i],
                   scale = scales[num.r],
                   AIC = results_road_all_n[i,num.r,1],
                   weight = akaike.weights(results_road_all_n[i,,1])$weights[num.r],
                   sex = unique(deer.i$sex),
                   age = unique(deer.i$age),
                   site = unique(deer.i$site),
                   season = ifelse(substr(unique(deer.i$id_yr_sn),13,15)=="win","win","sum"),
                   hr.area = area.n[i],
                   f.hr = f.vals.n[i],
                   r.hr = r.vals.n[i])
  scale.road.df <- rbind(scale.road.df,df)
}

## Roads South--------------#
deer.names <- unique(deer.s$id_yr_sn)
for (i in 1:length(deer.names)) {
  num.r <- which.min(results_road_all_s[i,,1])
  deer.i <- deer.s[deer.s$id_yr_sn==deer.names[i],]
  df <- data.frame(id = deer.names[i],
                   scale = scales[num.r],
                   AIC = results_road_all_s[i,num.r,1],
                   weight = akaike.weights(results_road_all_s[i,,1])$weights[num.r],
                   sex = unique(deer.i$sex),
                   age = unique(deer.i$age),
                   site = unique(deer.i$site),
                   season = ifelse(substr(unique(deer.i$id_yr_sn),13,15)=="win","win","sum"),
                   hr.area = area.s[i],
                   f.hr = f.vals.s[i],
                   r.hr = r.vals.s[i])
  scale.road.df <- rbind(scale.road.df,df)
}

## Forest North--------------#
scale.for.df <- data.frame()
deer.names <- unique(deer.n$id_yr_sn)
for (i in 1:length(deer.names)) {
  num.r <- which.min(results_for_n[i,,1])
  deer.i <- deer.n[deer.n$id_yr_sn==deer.names[i],]
  df <- data.frame(id = deer.names[i],
                   scale = scales[num.r],
                   AIC = results_for_n[i,num.r,1],
                   weight = akaike.weights(results_for_n[i,,1])$weights[num.r],
                   sex = unique(deer.i$sex),
                   age = unique(deer.i$age),
                   site = unique(deer.i$site),
                   season = ifelse(substr(unique(deer.i$id_yr_sn),13,15)=="win","win","sum"),
                   hr.area = area.n[i],
                   f.hr = f.vals.n[i],
                   r.hr = r.vals.n[i])
  scale.for.df <- rbind(scale.for.df,df)
}

## Forest South--------------#
deer.names <- unique(deer.s$id_yr_sn)
for (i in 1:length(deer.names)) {
  num.r <- which.min(results_for_s[i,,1])
  deer.i <- deer.s[deer.s$id_yr_sn==deer.names[i],]
  df <- data.frame(id = deer.names[i],
                   scale = scales[num.r],
                   AIC = results_for_s[i,num.r,1],
                   weight = akaike.weights(results_for_s[i,,1])$weights[num.r],
                   sex = unique(deer.i$sex),
                   age = unique(deer.i$age),
                   site = unique(deer.i$site),
                   season = ifelse(substr(unique(deer.i$id_yr_sn),13,15)=="win","win","sum"),
                   hr.area = area.s[i],
                   f.hr = f.vals.s[i],
                   r.hr = r.vals.s[i])
  scale.for.df <- rbind(scale.for.df,df)
}


scale.for.df <- read.csv("scale_forest_data.csv",header = T)
scale.road.df <- read.csv("scale_road_data.csv",header = T)

#--------------------------------------#
####   F O R E S T   M O D E L S    ####
#--------------------------------------#
# Global model
fmod.all <- lm(scale ~ site + scale(log(hr.area)) + scale(r.hr) + scale(f.hr) + 
                 sex + age + season, data = scale.for.df, weights = weight, na.action = "na.fail")
summary(fmod.all)


#-------------------------------------#
####     R O A D   M O D E L S     ####
#-------------------------------------#
# Global model
rmod.all <- lm(scale ~ sex + age + site + season + scale(log(hr.area)) + 
                 scale(r.hr) + scale(f.hr), data = scale.road.df, weights = weight, na.action = "na.fail")
summary(rmod.all)






