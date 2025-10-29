# Welcome to the ground.predator.viewshed() function download!
Build viewsheds from digital surface model and digital terrain model inputs in R to approximate the amount of cover in a given area.
To learn about potential research avenues check out: Doixx



![Untitled-1](https://github.com/user-attachments/assets/050c8372-c81a-4646-8f23-1ede740921d4)

Clockwise, A) pronghorn are a prey species in the Yellowstone ecosystem and spend time in wide expansive grasslands with a matrix of hills, and vegetation structures made of shrubs and trees which can restrict viewing the greater area. B) Coefficient graph of best fitting model explaining pronghorn movement and habitat selection after controlling for escape terrain and forage availability. C) Distribution of percent of visible area of used and available points in analysis split between bare earth and vegetation derived viewshed using custom function. Used and available values are relatively similar in both DTM derived distributions (red), but there is a distributional difference in the DSM derived values between the used and available locations. D) A representation of the benefits of using a custom function to calculate viewsheds of terrestrial predators using a digital surface and terrain model together. E) A representation of an output of percent cover area within a 500 meter buffer. As seen in the digital surface model, the viewshed accounts for vegetation that could hide an animal, unless the target height is greater than the vegetation structure. F) A DTM alone does not capture structural characteristics (top), a DSM alone does not capture structural characteristics that would mask a predator at a given height and over estimates visible areas where animals could not access such as the top of trees (middle), and finally a custom function output that accounts for a target height derived from the DTM model, but a line of sight calculation from the DSM. Pronghorn image: Wikipedia / Tobias Klenze / CC-BY-SA 4.0.

# How to implement with parallel and single core set up.

**Inputs:**

**buffersize** *Must be a numeric value. This is the radius distance used to create the circular buffer which the viewshed will be run on.*

**observer** *Must be two numeric values in a c(). These are the x,y coords of the observer location. Must match the crs of the dsm and dtm.*

**dsm** *Input equal area crs Digital Surface Model. Must be a terra SpatRaster object, must be able to be stacked with dtm ie. same resolution and crs(). Must also have a cell size that is a square for matrix math.*

**dtm** *Input equal area crs Digital Terrain Model. Must be a terra SpatRaster object, must be able to be stacked with dsm ie. same resolution and crs(). Must also have a cell size that is a square for matrix math.*

**obs_h** *Observer height, must be numeric. Unit is in the same as the elevation in the dsm/dtm, this is the value that is added to the p location to get z.*

**pred_h** *Predator height, must be numeric. Unit is in the same as the elevation in the dsm/dtm, this is the value that is added to each target location to draw line of sight between p's z location and target location + pred_h.*

**Outputs:** 

*Output of function is a terra SpatRaster with a binary (1,0) value. 1 is a visible sight based on parameters. 0 is not visible based on parameters.*

# On a single core
```r
#Running on a single core
#packages
require(terra)

#load custom function onto machine
custom_cpp <-  'https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/refs/heads/main/ground.predator.viewshed_forgithub.R'
source(custom_cpp)

# This function 1) Installs required RCpp package, 2) loads the custom C++ script into your temporary directory,
# and 3) compiles the C++ script for your machine. Internet is required to download. 
#This must be run every new session akin to a library()
ground.predator.viewshed.load()

#using terra, define input rasters...
r_dsm <- rast("/dsm.tif")
r_dtm <- rast("/dtm.tif")

#r_dsm and r_dtm must be stackable
r <- c(r_dsm,r_dtm)

#observer must be a c(x,y) in the coordinate system of the dsm and dtm
#x and y must be numeric.
p <- cbind(x,y)

#output is binary (0,1) raster, 1 is visible, 0 is not visible.
view_dsm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=dsm_yell,dtm=dsm_yell,obs_h=1.2,pred_h=1.0)

#optional, if you want to run as a bare earth viewshed
view_dtm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=r_dtm,dtm=r_dtm,obs_h=1.2,pred_h=1.0)

#Output is a terra SpatRaster
view_dsm

#Output as a polygon
as.poly(view_dsm)

```

# As a parallel process
```r
#Running in parallel
#packages
require(terra)
require(doParallel)
require(foreach)

#prepare workers
num_cores <- parallel::detectCores() - 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#compile for each worker the C++ code prior to running foreach loop
clusterEvalQ(cl, {
  
  # This function 1) Installs required RCpp package, 2) loads the custom C++ script into the worker temporary directory,
  # and 3) compiles the C++ script for your machine. Internet is required to download. 
  custom_cpp <-  'https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/refs/heads/main/ground.predator.viewshed_forgithub.R'
  source(custom_cpp)
  ground.predator.viewshed.load()
})

#remove any main system copies of function to not send to each worker. Must get unique compile in each worker...
rm("computeViewshedOptimized","ground.predator.viewshed","ground.predator.viewshed.load")

#run in parallel
foreach(.packages = c("terra", "sf","dplyr")) %dopar% {
  
  #using terra, define input rasters...
  r_dsm <- rast("/dsm.tif")
  r_dtm <- rast("/dtm.tif")
  
  #r_dsm and r_dtm must be stackable
  r <- c(r_dsm,r_dtm)
  
  #observer must be a c(x,y) in the coordinate system of the dsm and dtm
  #x and y must be numeric.
  p <- cbind(x,y)
  
  #output is binary (0,1) raster, 1 is visible, 0 is not visible.
  view_dsm <- ground.predator.viewshed(buffersize=500,observer=p,dsm=r_dsm,dtm=r_dtm,obs_h=1.2,pred_h=1.0)
  
  }

```

