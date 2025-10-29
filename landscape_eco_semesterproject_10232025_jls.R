#Entire workflow start to finish of Landscape Ecology semester project
#Jayden Skelly, MS student, Montana State University, jaydenskelly@montana.edu, Fall 2025
## ------------------------------------------------------------------------------------------
#This script companions my final project and presentation examining whether viewsheds can improve
#behavioral ecology movement modeling in a predator landscape using a case study of pronghorn in Yellowstone
## ------------------------------------------------------------------------------------------
#This work is sectioned out into:
#Spatial data prep
#Movement data prep, crop movement data within lidar area prep, availablilty prep
#Prepare random sample locations and draw availability for SSF
#Extraction of variables including calculating viewsheds
#Fitting SSF
#Compare AIC between terrain, VRM, ruggedness, viewshed-terrain, viewshed-surface while controlling for biomass/roughness
#Generate relevant output tables, graphs, maps
## ------------------------------------------------------------------------------------------
#Packages
#for spatial work
require(tidyverse)
require(sf)
require(terra)
require(mapview)
#for movement work
require(lubridate)
require(devtools)
#install_github("jmerkle1/MerkleLab-Code-Repository", subdir="MoveTools")
require(MoveTools)
require(doParallel)
require(foreach)
#for viewshed
custom_cpp <-  'https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/refs/heads/main/ground.predator.viewshed_forgithub.R'
source(custom_cpp)
ground.predator.viewshed.load()

#for analysis
require(MuMIn)
require(lme4)
require(MASS)
#require(caret) #masks cluster from survival, kept off unless needed
require(survival)
require(AICcmodavg)
#plotting
require(patchwork)
require(broom.mixed)
require(GGally)
require(ggrepel)
## ------------------------------------------------------------------------------------------
#Set working directory
setwd("/home/c32v191")

## ------------------------------------------------------------------------------------------
#Larger functions

## ------------------------------------------------------------------------------------------
#Mosaic NEON DSM and DTM into single tifs
dir_path <- "./Class/project/input_data/NEON_yell_2020/"
tif_files <- list.files(path = dir_path, pattern = "dsm.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
rasters <- lapply(tif_files, rast)
dsm_yell <- do.call(mosaic, rasters)

###
dir_path <- "./Class/project/input_data/NEON_yell_2020/"
tif_files <- list.files(path = dir_path, pattern = "dtm.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
rasters <- lapply(tif_files, rast)
dtm_yell <- do.call(mosaic, rasters)

plot(dsm_yell)
plot(dtm_yell)

res(dsm_yell)

writeRaster(dsm_yell, "./Class/project/output_data/dsm_yell.tif",overwrite=TRUE)
writeRaster(dtm_yell, "./Class/project/output_data/dtm_yell.tif",overwrite=TRUE)

#or 
dsm_yell <- rast("./Class/project/output_data/dsm_yell.tif")
dtm_yell <- rast("./Class/project/output_data/dtm_yell.tif")


#next create a buffer of the dsm/dtm, to use and clip out pronghorn points. 
#this will be determined by taking the outside edge, and subtracting into the raster the flight distance to be used
#in this case I will use a flight distance of 500 meters.

r_modified <- app(dtm_yell, fun=function(x) {
  x[!is.na(x)] <- 1
  return(x)
})
plot(r_modified)
r_modified <-  aggregate(r_modified, fact = 10, fun = mean)
res(r_modified)
bounds <- as.polygons(r_modified, dissolve = TRUE)
plot(bounds)
l <- as.lines(bounds)
buff_bounds <- buffer(l, 500, quadsegs=10, capstyle="round", 
       joinstyle="round", mitrelimit=NA, singlesided=FALSE)
plot(buff_bounds)

#clip bounds by the buff_bounds to get a smaller polygon in fact 500m smaller
bounds_steady <- crop(bounds,buff_bounds)
bounds_steady <- erase(bounds, bounds_steady)
plot(bounds_steady)

mapview(bounds)+mapview(l)+mapview(bounds_steady)

#save for later
saveRDS(bounds_steady, "./Class/project/output_data/bounds_steady.rds")
rm(buff_bounds,l,r_modified,bounds)   

#or
bounds_steady <- readRDS("./Class/project/output_data/bounds_steady.rds")
## ------------------------------------------------------------------------------------------
#Prepare Biomass
biomass_20 <- rast("./Class/project/input_data/rapr_biomass_2020_ouput.tif")
biomass_20 <- crop(biomass_20,project(bounds_steady,crs(biomass_20)))
writeRaster(biomass_20,"./Class/project/output_data/biomass_20.tif")

biomass_21 <- rast("./Class/project/input_data/rapr_biomass_2021_ouput.tif")
biomass_21 <- crop(biomass_21,project(bounds_steady,crs(biomass_21)))
writeRaster(biomass_21,"./Class/project/output_data/biomass_21.tif")

biomass_22 <- rast("./Class/project/input_data/rapr_biomass_2022_ouput.tif")
biomass_22 <- crop(biomass_22,project(bounds_steady,crs(biomass_22)))
writeRaster(biomass_22,"./Class/project/output_data/biomass_22.tif")

#or
biomass_20 <- rast("./Class/project/output_data/biomass_20.tif")
biomass_21 <- rast("./Class/project/output_data/biomass_21.tif")
biomass_22 <- rast("./Class/project/output_data/biomass_22.tif")

## ------------------------------------------------------------------------------------------
#Prepare tree cover
treecov <- rast('./Class/project/input_data/treecov.tif')
res(treecov)
treecov <- ifel(is.na(treecov), 0, treecov)
writeRaster(treecov,"./Class/project/output_data/treecov.tif")

#or
treecov <- rast("./Class/project/output_data/treecov.tif")

## ------------------------------------------------------------------------------------------
#Create VRM, Ruggedness, TPI
r <-  aggregate(dtm_yell, fact = 30, fun = mean)
res(r)

#TPI, TRI, TRIriley, TRIrmsd
rough <- terra::terrain(r, v = "roughness", neighbors = 8)
TPI <- terra::terrain(r, v = "TPI", neighbors = 8)
TRIriley <- terra::terrain(r, v = "TRIriley", neighbors = 8)
TRIrmsd <-terra::terrain(r, v = "TRIrmsd", neighbors = 8)
r <- c(r,rough,TPI,TRIriley,TRIrmsd)
names(r)[1] <- "DTM"
writeRaster(r,"./Class/project/output_data/r.tif")

#or
r <- rast("./Class/project/output_data/r.tif")
## ------------------------------------------------------------------------------------------
#SECTION 2: movement clean up, sample availability at each step, extract spatial data for movement and random SRS of landscape
## ------------------------------------------------------------------------------------------
#Import gps movement data, pre cleaned...
data <- readRDS("./Class/project/input_data/cleaned_pronghorndata_05292025.rds")
head(data)

#filter to Paradise
data <- data %>% filter(Herd == 'Paradise')

#filter down to points within the boundary
bounds_steady <- st_as_sf(bounds_steady)
bounds_steady <- st_transform(bounds_steady, crs = crs(data))
data <- st_filter(data, bounds_steady, .predicate = st_intersects)
mapview(data)

names(data)
unique(data$id_yr_seas)

#Quite a few individuals! I need to re build the step flags probably very intermitent
#I will plan to use the fine scale 1 hour data...
ids <- data %>% group_by(id_yr_seas) %>% summarise(n_rows = n()) %>% arrange(n_rows)
ids_keep <- ids %>% filter(n_rows > 100) %>% pull(id_yr_seas)
data <- data %>% filter(id_yr_seas %in% ids_keep)

#drop old movement columns
names(data)
data <- data %>% select(-c(abs.angle,dist,dt,speed,rel.angle,StepFlag))


data <- data %>% 
  arrange(id_yr_seas, DateTimeLocal)

#Calc Move Params
data$burst <- CalcBurst(data=data, id = TRUE, id_name="id_yr_seas", 
                         date_name="DateTimeLocal", Tmax = 3660*2) #set Tmax to a little more than double the fix rate, 3600 = 1 hour

data <- CalcMovParams(data=data, id_name = "id_yr_seas", 
                       date_name = "DateTimeLocal", burst=data$burst)
head(data)
saveRDS(data, "./Class/project/input_data/pronghorn_ready_landscapeeco.rds")

#or
data <- readRDS("./Class/project/input_data/pronghorn_ready_landscapeeco.rds")
## ------------------------------------------------------------------------------------------
#Then sample availability

# Take out of sf, and keep simply as a dataframe
proj <- st_crs(data)   # grab the projection for later
data <- data %>% 
  mutate(x = st_coordinates(.)[,1], # add x and y as columns
         y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry()  # need to remove the geometry column
head(data)
class(data)

hist(data$dt/3600)
table(data$dt)   # check to make sure all your steps are the same or very similar length
table(data$dt/3600) # now in hrs   
# Note, because of small variations in when the collar actually takes a location,
# you may see a variety of values here. All your steps should be relatively close in dt.
# For example, if you have 3 hour data, you only want to use steps that 
# are between 2.8 and 3.2 hours apart. 

# This is where you say 'I don't want these steps connected for further development of the SSF'
# Anytime there is a FALSE in the StepFlag column, the 'step sampling' code knows NOT to include the step
# you will need to change the next two lines of code !!!!
sum(data$StepFlag)   # this is how many 'actual/usable' steps you currently have
data$StepFlag[data$dt < 3600*0.75] <- FALSE   # this says 'don't connect steps less than 3/4 hrs apart
data$StepFlag[data$dt > 3600*1.25] <- FALSE   # this says 'don't connect steps greater than 1 1/4 hrs apart
sum(data$StepFlag)   # this is how many 'actual/usable' steps you now have

table(data$dt[data$StepFlag == TRUE])

numb.avail <- 15

#NOTE: Due to running this on the HPC a default paralell section of the function was causing issues
#to resolve this I created and stored the custom version which request soley 4 cores every time. It can be run at bottom of script

data_samp <- DrawRandSteps_emperical_custom(data=data, nr=numb.avail,   # nr = number of random steps 
                                     simult=FALSE,        # should it take a step length and a turning angle simultaneously?
                                     id_name="id_yr_seas", date_name="DateTimeLocal", x_name="x", y_name="y",   # what are the names of your id and date columns?
                                     withinID=TRUE,        # should sample from same individual or if FALSE all individuals in dataset?
                                     distMax = Inf, uniform.angles=FALSE)


data_samp <- drop_na(data_samp,"x_end")
data_samp <- st_as_sf(data_samp, coords=c("x_end","y_end"), dim="XY", crs=5072)
head(data)

#filter down to points within the boundary
bounds_steady <- st_as_sf(bounds_steady)
bounds_steady <- st_transform(bounds_steady, crs = crs(data_samp))
data_samp <- st_filter(data_samp, bounds_steady, .predicate = st_intersects)
#mapview(data_samp)

saveRDS(data_samp, "./Class/project/input_data/pronghorn_ready_landscapeeco_sampled.rds")

#or
data <- readRDS("./Class/project/input_data/pronghorn_ready_landscapeeco_sampled.rds")
## ------------------------------------------------------------------------------------------
#Now extract
data$DTM <- terra::extract(r$DTM, st_transform(data, crs=st_crs(r)), ID=FALSE)[,1]
data$roughness <- terra::extract(r$roughness, st_transform(data, crs=st_crs(r)), ID=FALSE)[,1]
data$TPI <- terra::extract(r$TPI, st_transform(data, crs=st_crs(r)), ID=FALSE)[,1]
data$TRIriley <- terra::extract(r$TRIriley, st_transform(data, crs=st_crs(r)), ID=FALSE)[,1]
data$TRIrmsd <- terra::extract(r$TRIrmsd, st_transform(data, crs=st_crs(r)), ID=FALSE)[,1]
data$treecov <- terra::extract(treecov, st_transform(data, crs=st_crs(treecov)), ID=FALSE)[,1]

data$BAFG[data$year_bio == 2020] <- terra::extract(biomass_20$annual_forb_and_grass, st_transform(data[data$year_bio == 2020, ], crs = st_crs(biomass_20)),ID = FALSE)[, 1]
data$BAFG[data$year_bio == 2021] <- terra::extract(biomass_21$annual_forb_and_grass, st_transform(data[data$year_bio == 2021, ], crs = st_crs(biomass_21)),ID = FALSE)[, 1]
data$BAFG[data$year_bio == 2022] <- terra::extract(biomass_22$annual_forb_and_grass, st_transform(data[data$year_bio == 2022, ], crs = st_crs(biomass_22)),ID = FALSE)[, 1]

data$BPFG[data$year_bio == 2020] <- terra::extract(biomass_20$perennial_forb_and_grass, st_transform(data[data$year_bio == 2020, ], crs = st_crs(biomass_20)),ID = FALSE)[, 1]
data$BPFG[data$year_bio == 2021] <- terra::extract(biomass_21$perennial_forb_and_grass, st_transform(data[data$year_bio == 2021, ], crs = st_crs(biomass_21)),ID = FALSE)[, 1]
data$BPFG[data$year_bio == 2022] <- terra::extract(biomass_22$perennial_forb_and_grass, st_transform(data[data$year_bio == 2022, ], crs = st_crs(biomass_22)),ID = FALSE)[, 1]

## ------------------------------------------------------------------------------------------
#Then calculate viewshed at each use/available point
#data <- data %>% head(600)


#transform to crs of rasters
crs = crs(dsm_yell)

data <- data %>% st_transform(crs)

#drop geometry
data <- data %>%
  mutate(
    x_end = st_coordinates(.)[, "X"],
    y_end = st_coordinates(.)[, "Y"]
  ) %>% st_drop_geometry()
                                              
num_cores <- parallel::detectCores() - (202)  #for 60cores, 300gb ram, -196 300/60 is 5gb a core #for 5 cores, 251
cl <- makeCluster(num_cores)
registerDoParallel(cl)


#compile for each worker the C++ code
clusterEvalQ(cl, {
  custom_cpp <-  'https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/refs/heads/main/ground.predator.viewshed_forgithub.R'
  source(custom_cpp)
  ground.predator.viewshed.load()
})

#remove any main system copies of function to not send...
rm("computeViewshedOptimized","ground.predator.viewshed","ground.predator.viewshed.load")
rm("biomass_20","biomass_21","biomass_22","bounds_steady","dsm_yell","dtm_yell","treecov","r")

#run in parallel
results <- foreach(i = seq_len(nrow(data)), .combine = rbind, .packages = c("terra", "sf","dplyr")) %dopar% {
  
  p <- cbind(data[i,]$x_end,data[i,]$y_end)

  dsm_yell <- rast("./Class/project/output_data/dsm_yell.tif")
  dtm_yell <- rast("./Class/project/output_data/dtm_yell.tif")
  
  view_dsm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dsm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)
  view_dtm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dtm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)

  eq_true_dsm <- sum(values(view_dsm), na.rm = TRUE)
  eq_total_dsm <- sum(!is.na(values(view_dsm)))
  p_view_dsm <- eq_true_dsm / eq_total_dsm

  eq_true_dtm <- sum(values(view_dtm), na.rm = TRUE)
  eq_total_dtm <- sum(!is.na(values(view_dtm)))
  p_view_dtm <- eq_true_dtm / eq_total_dtm
  
  c(p_view_dsm, p_view_dtm)
}

data$p_view_dsm <- results[, 1]
data$p_view_dtm <- results[, 2]
saveRDS(data, "./Class/project/input_data/pronghorn_ready_landscapeeco_extracted.rds")

#or
data <- readRDS("./Class/project/input_data/pronghorn_ready_landscapeeco_extracted.rds")
## ------------------------------------------------------------------------------------------
data_sensitivity <- data

#randomly select 50 strata to compare percent area view across different buffer size
ids <- data_sensitivity %>% filter(case==1) %>% sample_n(50) %>% pull(strata)
data_sensitivity <- data_sensitivity %>% filter(strata %in% ids)

#buffer sizes 50, 100, 150, 200, 300, 400...
buffersizes <- c(50,100,150,200,300,400)
for(buffer_s in buffersizes){
  print(buffer_s)


  num_cores <- parallel::detectCores() - (202)  #for 60cores, 300gb ram, -196 300/60 is 5gb a core #for 5 cores, 251
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  
  #compile for each worker the C++ code
  clusterEvalQ(cl, {
    custom_cpp <-  'https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/refs/heads/main/ground.predator.viewshed_forgithub.R'
    source(custom_cpp)
    ground.predator.viewshed.load()
  })
  
  #remove any main system copies of function to not send...
  rm("computeViewshedOptimized","ground.predator.viewshed","ground.predator.viewshed.load")
  rm("biomass_20","biomass_21","biomass_22","bounds_steady","dsm_yell","dtm_yell","treecov","r")
  
  #run in parallel
  results <- foreach(i = seq_len(nrow(data_sensitivity)), .combine = rbind, .packages = c("terra", "sf","dplyr")) %dopar% {
    
    p <- cbind(data_sensitivity[i,]$x_end,data_sensitivity[i,]$y_end)
    
    dsm_yell <- rast("./Class/project/output_data/dsm_yell.tif")
    dtm_yell <- rast("./Class/project/output_data/dtm_yell.tif")
    
    view_dsm <- ground.predator.viewshed(buffersize=buffer_s,observer=p,dsm=dsm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)
    view_dtm <- ground.predator.viewshed(buffersize=buffer_s,observer=p,dsm=dtm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)
    
    eq_true_dsm <- sum(values(view_dsm), na.rm = TRUE)
    eq_total_dsm <- sum(!is.na(values(view_dsm)))
    p_view_dsm <- eq_true_dsm / eq_total_dsm
    
    eq_true_dtm <- sum(values(view_dtm), na.rm = TRUE)
    eq_total_dtm <- sum(!is.na(values(view_dtm)))
    p_view_dtm <- eq_true_dtm / eq_total_dtm
    
    c(p_view_dsm, p_view_dtm)
  }
  saveRDS(results, paste0("./Class/project/input_data/results_",buffer_s,".rds"))

  } #end of sensitivity for loop

#then add the different results to data sensitivity as new columns
results_50 <- readRDS("./Class/project/input_data/results_50.rds")
results_100 <- readRDS("./Class/project/input_data/results_100.rds")
results_150 <- readRDS("./Class/project/input_data/results_150.rds")
results_200 <- readRDS("./Class/project/input_data/results_200.rds")
results_300 <- readRDS("./Class/project/input_data/results_300.rds")
results_400 <- readRDS("./Class/project/input_data/results_400.rds")

data_sensitivity$p_view_dsm_50 <- results_50[, 1]
data_sensitivity$p_view_dtm_50 <- results_50[, 2]

data_sensitivity$p_view_dsm_100 <- results_100[, 1]
data_sensitivity$p_view_dtm_100 <- results_100[, 2]

data_sensitivity$p_view_dsm_150 <- results_150[, 1]
data_sensitivity$p_view_dtm_150 <- results_150[, 2]

data_sensitivity$p_view_dsm_200 <- results_200[, 1]
data_sensitivity$p_view_dtm_200 <- results_200[, 2]

data_sensitivity$p_view_dsm_300 <- results_300[, 1]
data_sensitivity$p_view_dtm_300 <- results_300[, 2]

data_sensitivity$p_view_dsm_400 <- results_400[, 1]
data_sensitivity$p_view_dtm_400 <- results_400[, 2]

#visually compare how percent area viewed change across buffer size
saveRDS(data_sensitivity, "./Class/project/input_data/data_sensitvity.rds")

#or
data_sensitivity <-  readRDS("./Class/project/input_data/data_sensitvity.rds")

#w
## ------------------------------------------------------------------------------------------
#Then create SRS of lidar area to calculate viewsheds and compare landscape patterns to TPI etc accounting for cover
bounds_steady <- st_as_sf(bounds_steady)
bounds_steady <- st_transform(bounds_steady, crs = crs(dsm_yell))
random_points <- st_sample(
  x = bounds_steady,
  size = 1000,
  type = "random",
  exact = TRUE
) 
random_points <- do.call(rbind, random_points) %>% as.data.frame()
random_points$x <- random_points$V1
random_points$y <- random_points$V2
random_points$V1 <- NULL
random_points$V2 <- NULL
#create buffers
dsm_p_meta <- c()
start = Sys.time()
 for (i in 1:nrow(random_points)){
  print(i)
  p <- cbind(random_points[i,]$x,random_points[i,]$y)
  view_dsm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dsm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)
  eq_true_dsm <- sum(values(view_dsm), na.rm = TRUE)
  eq_total_dsm <- sum(!is.na(values(view_dsm)))
  p_view_dsm <- eq_true_dsm / eq_total_dsm
  dsm_p_meta <- rbind(dsm_p_meta, p_view_dsm)
}
random_points$p_view_dsm <- dsm_p_meta[, 1]
end = Sys.time()
total = start - end
random_points <- random_points %>% st_as_sf(coords=c(x='x',y='y'),crs=crs(dsm_yell))

#then extract the tree cover, TPI, ruggedness for all random samples...
r <- rast("./Class/project/output_data/r.tif")
random_points$DTM <- terra::extract(r$DTM, st_transform(random_points, crs=st_crs(r)), ID=FALSE)[,1]
random_points$roughness <- terra::extract(r$roughness, st_transform(random_points, crs=st_crs(r)), ID=FALSE)[,1]
random_points$TPI <- terra::extract(r$TPI, st_transform(random_points, crs=st_crs(r)), ID=FALSE)[,1]
random_points$TRIriley <- terra::extract(r$TRIriley, st_transform(random_points, crs=st_crs(r)), ID=FALSE)[,1]
random_points$TRIrmsd <- terra::extract(r$TRIrmsd, st_transform(random_points, crs=st_crs(r)), ID=FALSE)[,1]

treecov <- rast("./Class/project/output_data/treecov.tif")
random_points$treecov <- terra::extract(treecov, st_transform(random_points, crs=st_crs(treecov)), ID=FALSE)[,1]

saveRDS(random_points, "./Class/project/input_data/random_points_extracted.rds")
rm(random_points_buffers,view_dsm,view_dtm,p,dsm_crop,dtm_crop,dsm_p_meta,dsm_yell,dtm_yell,buffers,buffer,bounds_steady,cl,combined_sf,treecov,r,results)

#or
random_points <- readRDS("./Class/project/input_data/random_points_extracted.rds")
## ------------------------------------------------------------------------------------------
#SECTION 3: Analysis: LM of TPI to viewshed + cover, and model selection for SSF
## ------------------------------------------------------------------------------------------
#Models to compare SRS:

hist(random_points$p_view_dsm) #non-normal
hist(random_points$treecov)  #non-normal
hist(random_points$TPI)  #normal

# Fit a negative binomial GLM
glm_null <- glm.nb(p_view_dsm ~ 1, data = random_points)
glm_view_tpi <- glm.nb(p_view_dsm ~ TPI, data = random_points)
glm_view_tree <- glm.nb(p_view_dsm ~ treecov, data = random_points)
glm_view_rough <- glm.nb(p_view_dsm ~ roughness, data = random_points)
glm_view_tpi_tree <- glm.nb(p_view_dsm ~ TPI + treecov, data = random_points)
glm_view_rough_tree <- glm.nb(p_view_dsm ~ roughness + treecov, data = random_points)

cand.set.random <- list(glm_null,
                  glm_view_tpi,
                  glm_view_rough,
                  glm_view_tree,
                  glm_view_rough_tree,
                  glm_view_tpi_tree)

modnames <- c("Null", 
              "Terrain Position Index",
              "Roughness",
              "Tree Cover %",
              "Roughness + Tree Cover %",
              "TPI + Tree Cover %")

a.random <- aictab(cand.set = cand.set.random, modnames = modnames)
write.csv(a.random, "a_random.csv")
#lol with the top 2 models less than 2 AIC apart. It appears tree is very important...
#how strong?

#Validation
train_control <- trainControl(method = "cv", number = 5)
cv_model <- train(
  p_view_dsm ~ treecov,
  data = random_points,
  method = "glm.nb",
  trControl = train_control
)

validation <- cv_model$results
write.csv(validation, "validation.csv")
#pR2 = 0.37, RMSE = 0.099

summary(glm_view_tree)
summary(glm_view_tpi_tree)

## ------------------------------------------------------------------------------------------
#For now I will move onto the behavioral hypotheses ie the main part....

data$dist[data$dist == 0] <- 0.0000001
#drop geom
data <- data %>%
  sf::st_drop_geometry()

nrow(data)

#drop NAs
data$BPFG <- NULL
data <- drop_na(data)

drop <- data %>% group_by(strata) %>% summarise(usedsteps = sum(case)) %>% filter(usedsteps == 0) %>% pull(strata)
data <- data %>% filter(!strata %in% drop)
#========================================================
#Null model for use in AIC:
null_model <- clogit(case ~ dist + log(dist) +
                        strata(strata) + cluster(id_yr_seas), 
                      x=TRUE, y=TRUE,  
                      method = "efron", data=data)
#========================================================
variables <- c("TPI","roughness", "BAFG","p_view_dsm","p_view_dtm","treecov")
correl <- data %>% 
  dplyr::select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
correl <- correl %>% as.data.frame()
write.csv(correl, "correl.csv")
#Obviously cant put both viewsheds in same model
#========================================================
#Account for biomass
data$BAFG.st <- scale(data$BAFG)
BAFG_model <-  clogit(case ~ dist + log(dist) + BAFG.st +
                        strata(strata) + cluster(id_yr_seas), 
                      x=TRUE, y=TRUE,  
                      method = "efron", data=data)
summary(BAFG_model)
#And roughness
data$roughness.st <- scale(data$roughness)
roughness_model <-  clogit(case ~ dist + log(dist) + roughness.st +
                        strata(strata) + cluster(id_yr_seas), 
                      x=TRUE, y=TRUE,  
                      method = "efron", data=data)
#Both
roughness_biomass_model <-  clogit(case ~ dist + log(dist) + roughness.st + BAFG.st +
                             strata(strata) + cluster(id_yr_seas), 
                           x=TRUE, y=TRUE,  
                           method = "efron", data=data)
model.sel(null_model,
          BAFG_model,
          roughness_model,
          roughness_biomass_model,
          rank="AIC")
#========================================================
#Fit model accounting for different viewscapes
data$p_view_dsm.st <- scale(data$p_view_dsm)
data$p_view_dtm.st <- scale(data$p_view_dtm)
data$p_view_lgdsm.st <- scale(log(data$p_view_dsm))
data$p_view_lgdtm.st <- scale(log(data$p_view_dtm))

view_dsm_model <-  clogit(case ~ dist + log(dist) +  p_view_dsm.st +
                        strata(strata) + cluster(id_yr_seas), 
                      x=TRUE, y=TRUE,  
                      method = "efron", data=data)
view_dtm_model <-  clogit(case ~ dist + log(dist) +  p_view_dtm.st +
                            strata(strata) + cluster(id_yr_seas), 
                          x=TRUE, y=TRUE,  
                          method = "efron", data=data)
view_biomass_rough_dsm_model <-  clogit(case ~ dist + log(dist) + BAFG.st + roughness.st + p_view_dsm.st +
                            strata(strata) + cluster(id_yr_seas), 
                          x=TRUE, y=TRUE,  
                          method = "efron", data=data)
view_biomass_rough_dtm_model <-  clogit(case ~ dist + log(dist) + BAFG.st + roughness.st + p_view_dtm.st +
                            strata(strata) + cluster(id_yr_seas), 
                          x=TRUE, y=TRUE,  
                          method = "efron", data=data)
view_biomass_rough_lgdsm_model <-  clogit(case ~ dist + log(dist) + BAFG.st + roughness.st + p_view_lgdsm.st +
                                    strata(strata) + cluster(id_yr_seas), 
                                  x=TRUE, y=TRUE,  
                                  method = "efron", data=data)
view_biomass_rough_lgdtm_model <-  clogit(case ~ dist + log(dist) + BAFG.st + roughness.st + p_view_lgdtm.st +
                                      strata(strata) + cluster(id_yr_seas), 
                                    x=TRUE, y=TRUE,  
                                    method = "efron", data=data)

model.sel(roughness_biomass_model,
          view_dsm_model,
          view_dtm_model,
          view_biomass_rough_dsm_model,
          view_biomass_rough_dtm_model,
          view_biomass_rough_lgdsm_model,
          view_biomass_rough_lgdtm_model,
          rank="AIC")
#Dsm is better than terrain alone, and terrain viewscape is barely better than biomass alone...
#Now account for roughness and biomass to better account for biomass and available escape terrain to confirm
#========================================================
data$treecov.st <- scale(data$treecov)
treecov_model <-  clogit(case ~ dist + log(dist) + BAFG.st + roughness.st + treecov.st +
                                            strata(strata) + cluster(id_yr_seas), 
                                          x=TRUE, y=TRUE,  
                                          method = "efron", data=data)
model.sel(roughness_biomass_model,
          view_dsm_model,
          view_dtm_model,
          view_biomass_rough_dsm_model,
          view_biomass_rough_dtm_model,
          view_biomass_rough_lgdsm_model,
          view_biomass_rough_lgdtm_model,
          treecov_model,
          rank="AIC")
cand.set.test <- list(null_model,
                      BAFG_model,
                      roughness_model,
                      roughness_biomass_model,
                      view_dsm_model,
                      view_dtm_model,
                      view_biomass_rough_dsm_model,
                      view_biomass_rough_dtm_model,
                      view_biomass_rough_lgdsm_model,
                      view_biomass_rough_lgdtm_model,          
                      treecov_model)

modnames <- c("Intercept only",
              "Biomass",
              "Roughness",
              "Controls (Biomass + Roughness)", 
              "Viewshed (DSM)",
              "Viewshed (DTM)",
              "Controls + Viewshed (DSM)",
              "Controls + Viewshed (DTM)",
              "Controls + log (Viewshed (DSM))",
              "Controls + log (Viewshed (DTM))",
              "Controls + Tree Cover % (at 30 M)")

a.test <- aictab(cand.set = cand.set.test, modnames = modnames)
write.csv(a.test, "a_test.csv")
#========================================================
#========================================================
#Appendix analysis using the 1,000 sample
data_sensitivity$dist[data_sensitivity$dist == 0] <- 0.0000001
data_sensitivity <- data_sensitivity %>%
  sf::st_drop_geometry()
nrow(data_sensitivity)
data_sensitivity$BPFG <- NULL
data_sensitivity <- drop_na(data_sensitivity)
drop <- data_sensitivity %>% group_by(strata) %>% summarise(usedsteps = sum(case)) %>% filter(usedsteps == 0) %>% pull(strata)
data_sensitivity <- data_sensitivity %>% filter(!strata %in% drop)

data_sensitivity$p_view_dsm.st <- scale(data_sensitivity$p_view_dsm)
data_sensitivity$biomass.st <- scale(data_sensitivity$BAFG)

data_sensitivity$p_view_dsm_50.st <- scale(data_sensitivity$p_view_dsm_50)
data_sensitivity$p_view_dsm_100.st <- scale(data_sensitivity$p_view_dsm_100)
data_sensitivity$p_view_dsm_150.st <- scale(data_sensitivity$p_view_dsm_150)
data_sensitivity$p_view_dsm_200.st <- scale(data_sensitivity$p_view_dsm_200)
data_sensitivity$p_view_dsm_300.st <- scale(data_sensitivity$p_view_dsm_300)
data_sensitivity$p_view_dsm_400.st <- scale(data_sensitivity$p_view_dsm_400)

view_dsm_model_50 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_50.st + biomass.st +
                              strata(strata) + cluster(id_yr_seas), 
                              x=TRUE, y=TRUE,  
                              method = "efron", data=data_sensitivity)
view_dsm_model_100 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_100.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
view_dsm_model_150 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_150.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
view_dsm_model_200 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_200.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
view_dsm_model_300 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_300.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
view_dsm_model_400 <-  clogit(case ~ dist + log(dist) +  p_view_dsm_400.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
view_dsm_model_500 <-  clogit(case ~ dist + log(dist) +  p_view_dsm.st + biomass.st +
                               strata(strata) + cluster(id_yr_seas), 
                             x=TRUE, y=TRUE,  
                             method = "efron", data=data_sensitivity)
aic <- model.sel(view_dsm_model_50,
                 view_dsm_model_100,
                 view_dsm_model_150,
                 view_dsm_model_200,
                 view_dsm_model_300,
                 view_dsm_model_400,
                 view_dsm_model_500,
                 rank="AIC")
aic


cand.set <- list(view_dsm_model_50,
                 view_dsm_model_100,
                 view_dsm_model_150,
                 view_dsm_model_200,
                 view_dsm_model_300,
                 view_dsm_model_400,
                 view_dsm_model_500)
modnames <- c("050 meter + biomass", 
             "100 meter + biomass",
             "150 meter + biomass",
             "200 meter + biomass",
             "300 meter + biomass",
             "400 meter + biomass",
             "500 meter + biomass")

a <- aictab(cand.set = cand.set, modnames = modnames)
write.csv(a, "a.csv")
#========================================================
## ------------------------------------------------------------------------------------------
#SECTION 4: Interpretation: Create graphs, maps, tables
## ------------------------------------------------------------------------------------------
#moving on!
#I'd like a graph of a location and the viewshed from each output
#And AIC tables
#Pairwise table of variables


new_data <- data.frame(
  dist = mean(data$dist),
  BAFG.st = mean(data$BAFG.st),
  p_view_lgdsm.st = seq(min(data$p_view_lgdsm.st), max(data$p_view_lgdsm.st), length.out = 100),
  strata = data$strata[1], 
  id_yr_seas = data$id_yr_seas[1]
)

pred <- predict(view_biomass_lgdsm_model, newdata = new_data, type = "risk")
plot(new_data$p_view_lgdsm.st, pred, type = "l", xlab = "Viewshed", ylab = "Predicted Probability")
#========================================================
#Get some summary stats for paper
used <-  data %>% filter(case==1)
avail <-  data %>% filter(case==0)
range(used$p_view_dsm)
range(avail$p_view_dsm)
range(used$p_view_dtm)
range(avail$p_view_dtm)


g1 <- ggplot(used) +
  geom_density(aes(x = p_view_dtm, fill = "DTM"), color = "black", alpha = 0.5) +
  geom_density(aes(x = p_view_dsm, fill = "DSM"), color = "black", alpha = 0.5) +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("DTM" = "#CB4154", "DSM" = "tan")) +
  theme(legend.position = "none") +
  ggtitle("Used Locations")
g2 <- ggplot(avail) +
  geom_density(aes(x = p_view_dtm, fill = "DTM"), color = "black", alpha = 0.5) +
  geom_density(aes(x = p_view_dsm, fill = "DSM"), color = "black", alpha = 0.5) +
  labs(x = "% of 500 Meter Buffer Visible", y = "") +
  scale_fill_manual(
    values = c("DTM" = "#CB4154", "DSM" = "tan"),
    name = "Viewshed Type", 
    labels = c("DTM" = "Viewshed % DTM", "DSM" = "Viewshed % DSM") 
  ) +
  ggtitle("Available locations")+
  theme(
    legend.position = c(0.8, 0.8), 
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.7)))
g1 / g2
#========================================================
#Now plot these with color coding...
pos <- position_jitterdodge(
  jitter.width = 0, 
  dodge.width = 0.1  
)

# Define custom colors, making sure the names match the `term` values
custom_colors <- c(
  "BAFG.st" = "cadetblue",
  "dist" = "darkorange3",
  "log(dist)" = "goldenrod",
  "p_view_lgdsm.st" = "tomato1",
  "roughness.st" = "slateblue"
)

# Define the custom order of legend items using the `term` names
custom_order <- c("roughness.st", "p_view_lgdsm.st", "log(dist)", "dist", "BAFG.st")

ggplot(coeffs, aes(x = estimate, y = term, color = term)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, position = pos) +
  geom_point(position = pos, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = custom_colors,
    breaks = custom_order, # Use the `breaks` argument to set the order
    labels = c(
      "roughness.st" = "Roughness scaled",
      "p_view_lgdsm.st" = "Total area viewed (DSM) scaled",
      "log(dist)" = "Log (Step length) control",
      "dist" = "Step length control", 
      "BAFG.st" = "Biomass scaled"
    )
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    color = "Term"
  ) + 
  theme(
    legend.position = c(0.9, 0.6), 
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50"),
    aspect.ratio = 1,
    axis.text.y = element_blank(),
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.7)))
#========================================================
#create .kmz of a viewshed for a single point:
k <- data.frame(y=44.93239170852187,x=-110.61359913397023) %>% st_as_sf(coords=c(x='x',y='y'),crs=4326) %>% 
                   st_transform(crs=crs(dsm_yell)) %>% 
                   mutate( x = st_coordinates(.)[, "X"],
                           y = st_coordinates(.)[, "Y"]) %>%
                   st_drop_geometry()
p <- cbind(k[1,]$x,k[1,]$y)
view_dsm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dsm_yell,dtm=dsm_yell,obs_h=1.2,pred_h=1.0)
view_dtm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dtm_yell,dtm=dtm_yell,obs_h=1.2,pred_h=1.0)

as.poly <- as.polygons(view_dsm)
as.poly <- as.poly[as.poly$lyr.1 == 1, ]
output_file <- "example_dsm_b.kml"
writeVector(as.poly, output_file, overwrite=TRUE, filetype="KML")

crop(dsm_yell)
plot(as.poly)
#========================================================
token <- "7552632568:AAG_uHdDybt8SPx9ihaV_1y4j8-WOt-z4mg"
chat_id <- "7726819829"
message <- "Landscape Ecology Project Completed! "
url <- paste0("https://api.telegram.org/bot", token, "/sendMessage?chat_id=", chat_id, "&text=", URLencode(message))
httr::GET(url)
#========================================================
#Appendix
DrawRandSteps_emperical_custom <- function (data = data, nr = 10, distMax = Inf, simult = FALSE, 
          id_name = "id", date_name = "date", x_name = "x", y_name = "y", 
          withinID = FALSE, uniform.angles = FALSE) 
{
  if (all(c("circular", "plyr", "parallel") %in% utils::installed.packages()[, 
                                                                             1]) == FALSE) 
    stop("You must install the following packages: circular, plyr, and snowfall")
  require(circular)
  require(plyr)
  require(parallel)
  if (inherits(data, "sf") == TRUE) 
    stop("data should be a data.frame")
  if (inherits(data, "data.frame") != TRUE) 
    stop("data should be a data.frame")
  if (any(colnames(data) == date_name) == FALSE) 
    stop(print("Your date_name is incorrect."))
  if (any(colnames(data) == id_name) == FALSE) 
    stop(print("Your id_name is incorrect."))
  if (any(colnames(data) == x_name) == FALSE) 
    stop(print("Your x_name is incorrect."))
  if (any(colnames(data) == y_name) == FALSE) 
    stop(print("Your y_name is incorrect."))
  if ("x_end" %in% names(data) | "strata" %in% names(data)) 
    stop("You have already calculated Random Steps! Need to start over with a clean dataframe to recalculate.")
  if (inherits(data[, date_name], "POSIXct") != TRUE) 
    stop("date_name is not POSIXct")
  if (any(duplicated(data[c(id_name, date_name)])) == TRUE) 
    stop("You have duplicates in your database")
  if (any(is.na(data[, date_name]) == TRUE)) 
    stop("You have NAs in your date column")
  if (any(is.na(data[, id_name]) == TRUE)) 
    stop("You have NAs in your id column")
  step_dt <- unique(na.omit(data$dt[data$StepFlag == TRUE]))
  if (diff(range(step_dt)) > 1200) 
    print("WARNING!!!!! You have a mix of step lengths in this dataset, and they differ by more than 20 mins!!! It is not a good idea to draw random steps without fixing this.")
  if (any(c("rel.angle", "abs.angle") %in% colnames(data) == 
          FALSE) == TRUE) 
    stop("you do not have an abs.angle and rel.angle columns in data. You need to rerun calcMovParams")
  rm(step_dt)
  key <- 1:nrow(data)
  key2 <- key[order(data[, id_name], data[, date_name])]
  if (all(key == key2) == FALSE) 
    stop(print("Your data are not ordered correctly"))
  rm(key, key2)
  u <- unique(data[, id_name])
  extr_col_nms <- names(data)
  extr_col_nms <- extr_col_nms[extr_col_nms %in% c(id_name, 
                                                   date_name, x_name, y_name, "burst", "dist", "dt", "speed", 
                                                   "abs.angle", "rel.angle", "StepFlag") == FALSE]
  data$strata <- 1:nrow(data)
  data$strata[data$StepFlag == FALSE] <- NA
  data$strata <- as.numeric(as.factor(data$strata))
  data$case <- ifelse(is.na(data$strata) == TRUE, NA, 1)
  data$x_end <- c(data[2:nrow(data), x_name], NA)
  data$y_end <- c(data[2:nrow(data), y_name], NA)
  data$x_end[is.na(data$strata) == TRUE] <- NA
  data$y_end[is.na(data$strata) == TRUE] <- NA
  north <- c(NA, data$abs.angle[1:(nrow(data) - 1)])
  no_cores <- 4
  clust <- parallel::makeCluster(no_cores)
  parallel::clusterExport(clust, varlist = c("data", "simult", 
                                             "withinID", "distMax", "uniform.angles", "nr", "u", "north", 
                                             "id_name", "date_name", "x_name", "y_name"), envir = environment())
  temp <- parallel::clusterApplyLB(clust, 1:length(u), function(i) {
    library(circular)
    num <- length(data$strata[is.na(data$strata) == FALSE & 
                                data[, id_name] == u[i]])
    if (num != 0) {
      if (withinID == FALSE) {
        dlist2 <- data[data$StepFlag == TRUE, ]
        dlist <- NA
      }
      else {
        dlist <- data[data$StepFlag == TRUE & data[, 
                                                   id_name] == u[i], ]
        dlist2 <- NA
      }
      if (simult == TRUE) {
        if (withinID == FALSE) {
          ts <- dlist2[sample(1:nrow(dlist2), nr * num, 
                              replace = TRUE), ]
        }
        else {
          ts <- dlist[sample(1:nrow(dlist), nr * num, 
                             replace = TRUE), ]
        }
      }
      else {
        if (withinID == FALSE) {
          ts <- data.frame(rel.angle = sample(dlist2[, 
                                                     "rel.angle"], nr * num, replace = TRUE), 
                           dist = sample(dlist2[, "dist"], nr * num, 
                                         replace = TRUE))
        }
        else {
          ts <- data.frame(rel.angle = sample(dlist[, 
                                                    "rel.angle"], nr * num, replace = TRUE), 
                           dist = sample(dlist[, "dist"], nr * num, 
                                         replace = TRUE))
        }
      }
      rng <- range(data$strata[data[, id_name] == u[i]], 
                   na.rm = T)
      ts$strata <- rep(rng[1]:rng[2], each = nr)
      ts$case <- 0
      north <- north[is.na(data$strata) == FALSE & data[, 
                                                        id_name] == u[i]]
      if (length(north) != num) 
        stop("you have problem 2")
      north <- rep(north, each = nr)
      angles <- ts$rel.angle
      ts$rel.angle <- ifelse(ts$rel.angle + north < 360, 
                             ts$rel.angle + north, (ts$rel.angle + north) - 
                               360)
      if (distMax != Inf) {
        ts$dist[ts$dist > distMax] <- distMax
      }
      if (uniform.angles == TRUE) {
        ts$rel.angle <- runif(length(ts$rel.angle), 0, 
                              359.99999999)
        angles <- ts$rel.angle
      }
      d <- data[is.na(data$strata) == FALSE & data[, id_name] == 
                  u[i], ]
      ts <- merge(ts[, c("strata", "case", "rel.angle", 
                         "dist")], d[, c("x", "y", date_name, "strata")])
      ts$id <- u[i]
      ts$x_end <- ts[, x_name] + (cos(circular::rad(ts$rel.angle)) * 
                                    ts$dist)
      ts$y_end <- ts[, y_name] + (sin(circular::rad(ts$rel.angle)) * 
                                    ts$dist)
      ts$rel.angle <- angles
      names(ts)[names(ts) == "id"] <- id_name
      return(ts)
    }
    else {
      return(NULL)
    }
  })
  parallel::stopCluster(clust)
  temp <- do.call(rbind, temp)
  temp <- plyr::rbind.fill(data, temp)
  temp <- temp[order(temp[, id_name], temp[, date_name], temp$case, 
                     decreasing = F), ]
  temp <- temp[is.na(temp$strata) == FALSE, ]
  temp2 <- temp[temp$case == 1, c("strata", "dt", "StepFlag", 
                                  extr_col_nms)]
  temp$dt <- NULL
  temp$StepFlag <- NULL
  if (length(extr_col_nms) > 0) {
    for (i in 1:length(extr_col_nms)) {
      temp[, extr_col_nms[i]] <- NULL
    }
  }
  temp <- merge(temp, temp2, all.x = TRUE)
  temp <- temp[order(temp[, id_name], temp[, date_name], temp$case, 
                     decreasing = F), ]
  temp$speed <- temp$dist/temp$dt
  abs.angle <- atan2((temp$y_end - temp[, y_name]), (temp$x_end - 
                                                       temp[, x_name]))
  temp$abs.angle <- as.numeric(circular::conversion.circular(circular::circular(abs.angle, 
                                                                                units = "radians"), units = "degrees", zero = pi/2, rotation = "clock", 
                                                             modulo = "2pi"))
  row.names(temp) <- 1:nrow(temp)
  temp$burst <- NULL
  temp$StepFlag <- NULL
  temp <- temp[c("case", setdiff(names(temp), "case"))]
  temp <- temp[c("strata", setdiff(names(temp), "strata"))]
  return(temp)
}
