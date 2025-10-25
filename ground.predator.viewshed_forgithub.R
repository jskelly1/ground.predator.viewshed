#Jayden Skelly, 10-25-2025
#jaydenskelly@montana.edu

#Beta version of calculating prey view sheds to terrestrial predators.

#Example: ground.predator.viewshed.load()
ground.predator.viewshed.load <- function(){
  required_packages <- c("Rcpp", "terra", "dplyr")
  load_or_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  invisible(lapply(required_packages, load_or_install))
  print("Downloading compile from Github: https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/4089e8ef6cb1e5a79ac1f1104fdf3e4e7c701abb/viewshed2.cpp")
  
  
  # Download C++ file to temp directory
  temp_cpp <- tempfile(fileext = ".cpp")
  cpp_url <- "https://raw.githubusercontent.com/jskelly1/ground.predator.viewshed/4089e8ef6cb1e5a79ac1f1104fdf3e4e7c701abb/viewshed2.cpp"
  
  download.file(cpp_url, destfile = temp_cpp, mode = "wb")
  
  # Compile the C++ code
  Rcpp::sourceCpp(temp_cpp)
  
  
  print("Success: C++ code compiled and packages loaded.") }


#I need inputs:
#buffersize = a buffer around the observer location in meters
#observer = c(x,y), crs is assumed to be the same as dsm
#dtm = dtm to use when computing the ground height of potential predator
#dsm = dsm to use when computing the vegetation height masking potential terrestrial predator
#obs_h = observer height above dtm at point(effectively z value, will be added to dtm), think of as eye level.
#pred_h = predator height above dtm at at a given point(effectively z value, will be added to dtm), think of as eye level of predator.

#output will be a binary (1,0) raster
#Example: viewshed <- ground.predator.viewshed(buffersize=50,observer=c(531547, 4974948),dsm=dsm,dtm=dtm,obs_h=1.2,pred_h=1.5)
ground.predator.viewshed <- function(buffersize,observer,dsm,dtm,obs_h,pred_h){

  # Check: buffer size must be numeric
  if (!is.numeric(buffersize)) {
    stop("Error: buffer size must be numeric.")
  }
  
  # Check: observer must be numeric and length 2
  if (!is.numeric(observer) || length(observer) != 2) {
    stop("Error: observer must be numeric and contain two values (e.g., c(x, y)).")
  }
  
  # Check: dsm must be a terra::SpatRaster
  if (!inherits(dsm, "SpatRaster")) {
    stop("Error: dsm must be a terra::SpatRaster.")
  }
  
  # Check: dtm must be a terra::SpatRaster
  if (!inherits(dtm, "SpatRaster")) {
    stop("Error: dtm must be a terra::SpatRaster.")
  }
  
  # Check: obs_h must be numeric
  if (!is.numeric(obs_h)) {
    stop("Error: obs_h must be numeric.")
  }
  
  # Check: pred_h must be numeric
  if (!is.numeric(pred_h)) {
    stop("Error: pred_h must be numeric.")
  }
  
  # Check: ground.predator.viewshed.load function is loaded
  if (!exists("ground.predator.viewshed.load")) {
    stop("Error: ground.predator.viewshed.load function is not loaded.")
  }
  
  # Check: dtm and dsm must be stackable
  if (!terra::compareGeom(dtm, dsm, stopOnError = FALSE)) {
    stop("Error: dtm and dsm do not match. Must be able to be stacked.")
  }
  
  # Check: CRS must match
  if (!terra::crs(dtm) == terra::crs(dsm)) {
    stop("Error: CRS of dtm and dsm must match.")
  }
  
  # Check: observer is within raster extent
  ext <- terra::ext(dsm)
  if (!(observer[1] >= ext$xmin && observer[1] <= ext$xmax &&
        observer[2] >= ext$ymin && observer[2] <= ext$ymax)) {
    stop("Error: Observation point is off the raster. Are the coordinates from the same CRS?")
  }
  
  
  # Check: projection is equal-area
  if (!grepl("equal", terra::crs(dsm), ignore.case = TRUE)) {
    warning("Warning: DSM or DTM may not be in an equal-area projection.")
  }
  
  print("Computing Viewshed")
  obs_coords <- matrix(observer, ncol = 2)
  observer_point <- vect(obs_coords, crs = crs(dtm))
  buffer_poly <- buffer(observer_point, width = buffersize)
  
  dtm_cropped <- crop(dtm, buffer_poly)
  dsm_cropped <- crop(dsm, buffer_poly)
  
  dtm_mat <- as.matrix(dtm_cropped)
  dsm_mat <- as.matrix(dsm_cropped)
  
  row_col <- sqrt(length(dtm_mat))
  middle <- row_col/2
  
  dtm_matrix <- matrix(dtm_mat, nrow = row_col, ncol = row_col, byrow = TRUE)
  dsm_matrix <- matrix(dsm_mat, nrow = row_col, ncol = row_col, byrow = TRUE)
  
  cx = middle 
  cy = middle 
  
  viewshed <-  computeViewshedOptimized(dsm_matrix,dtm_matrix, cy, cx,obs_height = obs_h, tgt_height = pred_h)

  ex <- ext(dtm_cropped)
  cr <- crs(dtm_cropped)
  r <- rast(viewshed, extent = ex, crs = cr) %>% mask(buffer_poly)
  print("Success.")
  return(r)
}