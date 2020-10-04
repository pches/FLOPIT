###################################################
# file: FLOPIT_function.R
#
###################################################
# Author and copyright: K. Joel Roop-Eckart, 2018
# The Pennsylvania State University
# kjr30@psu.edu
#
# Distributed under the GNU general public license
# No warranty
#
###################################################
# Last changes: November, 2018 (K. Joel Roop-Eckart)
# Estimated run time on a single core: 30 minutes
#
###################################################
# load packages
require('raster')
require('rgdal')
require('sp')
require('FNN')
require('dismo')
require('deldir')
require('rgeos')
require('RColorBrewer')
require('scales')
require('gstat')

################################### FLOPIT function ##########################################
FLOPIT <- function(flood_rasters_names, flood_probabilities, elevation_raster, depth, aggregation_value, method, map_type,
                    save, save_path_data, save_path_map, save_path_zones){
################ Import Houston Clipped and Resampled Raster files ###########################
# be sure to change file paths
n_rp <- length(flood_rasters_names)
flood_wse_vector <- vector(length = n_rp)
# calculate return period values for each flood
return_periods <- 1/flood_probabilities # 4% flood and mean wse were removed due to data quality
probs_pct <- 100 * (1/return_periods) # probabilities associated with the return periods

# read flood rasters 
for(i in 1:n_rp){
    flood_wse_vector[i] <- list(raster(flood_rasters_names[i]))
    names(flood_wse_vector[[i]])<-'rasterdata'
}

# import Lidar DEM of the study area
dem <- raster(elevation_raster)
print('Done reading flood rasters and elevation data')

######################### Convert flood depths (ft) to elevations (ft)
if(depth==TRUE){
    # convert depths to water surface elevations
    for (i in 1:n_rp) {
        tmp_rast_vals <- getValues(flood_wse_vector[[i]])+getValues(dem)
        flood_wse_vector[i] <- list(setValues(flood_wse_vector[[i]], tmp_rast_vals))
    }
}
print('Done with calculating water surface elevation')

######################### Interpolate missing values
# aggregate the rasters
returnLevels_wse_upscaled <- vector(length = n_rp)

for(i in 1:n_rp){
  returnLevels_wse_upscaled[i] <- list(aggregate(flood_wse_vector[[i]], fact = aggregation_value, fun = mean))
}

dem_upscaled <- aggregate(dem, 
                      fact = aggregation_value, 
                      fun = mean)
print('Done with aggregating elevation and water surface elevation data')

######################### Smooth water surface elevation data using Inverse Weighted Distance method
# define inverse weighted distance interpolation/extrapolation function
# define raster grid that covers entire interpolation/extrapolation space (DEM suggested for this use)
flood_raster_grid <- rasterToPoints(dem_upscaled, spatial = TRUE)
gridded(flood_raster_grid) = TRUE

# define function
IWD_interp <- function(raster,grid){
    raster_points <- rasterToPoints(raster, spatial = TRUE)
    raster_gstat <- gstat(formula = rasterdata ~ 1, 
                          data = raster_points, 
                          nmax = 10, # for local kriging: the number of nearest observations that should be used for a kriging prediction or simulation
                          set = list(idp = 0.5)
                          )
    z <- predict(rast_gstat, grid)
    raster@data@values <- z@data[,1]
    return(raster)
  }

# start spatial WSE interpolation/extrapolation
returnLevels_wse_smooth <- vector(length = n_rp)
for (i in 1:n_rp){
  returnLevels_wse_smooth[i] <- list(IWD_interp(raster = returnLevels_wse_upscaled[[i]], 
                                                grid = flood_raster_grid)
                                     )}

# Test for NA values after interpolation: If there are any, interpolation has made a mistake
for (i in 1:n_rp) {
  out <- which(is.na(returnLevels_wse_smooth[[i]]@data@values)==TRUE)
  if(length(out)>0)(print('ERROR: Water surface elevation interpolation failed'))
}
print('Done with smooting water surface elevation')

######################### Interpolate flood probability map
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map <- dem_upscaled
# set all values of the new raster to be NA
setValues(rt_map, NA)
# define a vector vals to write interpolated flood return periods to
vals <- vector(mode = 'numeric', length = length(rt_map@data@values))
# for loop uses splines to define flood elevation to return period for each cell
# and interpolate the return period associated with each cell's elevation
elevations <- getValues(dem_upscaled)
pb <- txtProgressBar(min = 0, max = length(getValues(rt_map)), initial = 0, char = '=', style = 1)
start <- Sys.time()
floods <- vector(length = n_rp)
for (i in 1:length(rt_map@data@values)) {
  setTxtProgressBar(pb, i)
  for (j in 1:n_rp) {floods[j] <- returnLevels_wse_smooth[[j]]@data@values[i]}
  floods <- sort(floods, decreasing = FALSE) # due to occasional error introduced by aggregation, 
  # on steep slopes, smaller floods may produce higher WSE elevations for a single tile due to aggregation.
  # While these errors reduce accuracy, sorting the floods allows the interpolation scheme to make a
  # realistic approximation or the relationship and finish the job
  if(method == 'log-linear')(vals[i] <- 1/((10^((
    (elevations[i]-max(floods[which(floods<elevations[i])]))*
      
      (log10(max(probs_pct[which(floods>elevations[i])]))-log10(min(probs_pct[which(floods<elevations[i])])))/
      
      (min(floods[which(floods>elevations[i])])-max(floods[which(floods<elevations[i])]))
    
  )+log10(min(probs_pct[which(floods<elevations[i])]))))/100)) # official FEMA log-linear interpolation formula
  else if(method == 'spline')(vals[i] <- spline(floods,return_periods,xout = elevations[i], method = 'hyman')$y)
}
close(pb)
end <- Sys.time()
end-start

######################### Post process interpolated flood probability data
# return periods cannot be negative, replace negative return period with the minimum return period in the range
vals<-replace(vals,which(vals<min(return_periods)),min(return_periods))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals<-replace(vals,which(vals>max(return_periods)),NA)
# if the wse of the largest flood is less than the elevation, set the pixel value to NA.
vals<-replace(vals, which(getValues(returnLevels_wse_smooth[[n_rp]])<getValues(dem_upscaled)), NA)
# if the wse of the a flood is greater than the elevation, and the return period is NA, 
# set the pixel value to the return period of that flood.

for (i in 1:n_rp) {
  vals<-replace(vals, intersect(which(getValues(returnLevels_wse_smooth[[i]])>getValues(dem_upscaled)), which(is.na(vals)==TRUE)), return_periods[i])
}

for (i in 1:n_rp) {
  vals<-replace(vals, intersect(which(is.na(getValues(returnLevels_wse_upscaled[[i]]))==FALSE), which(vals>return_periods[i])), return_periods[i])
}

vals<-replace(vals, intersect(which(is.na(getValues(returnLevels_wse_upscaled[[n_rp]]))==FALSE), which(is.na(vals)==TRUE)), return_periods[n_rp])

vals<-replace(vals, which(is.na(getValues(returnLevels_wse_upscaled[[n_rp]]))==TRUE), NA)
# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in 
# the 500 year zone in the 100 year zone

for (i in 1:(n_rp-1)) {
  vals<-replace(vals, intersect(intersect(which(is.na(getValues(returnLevels_wse_upscaled[[i+1]]))==FALSE), 
                                          which(is.na(getValues(returnLevels_wse_upscaled[[i]]))==TRUE)),
                                which((vals<return_periods[i]))), return_periods[i])
}


#which(is.na(getValues(returnLevels_wse_upscaled[[i]]))==TRUE)

########################## recreate flood zones
# create new raster
floodzones <- dem_upscaled
# create vector where all land above 500 yr flood is NA
bounds500 <- getValues(returnLevels_wse_upscaled[[which(return_periods==500)]])
# replace all non NA values with 500
bounds500 <- replace(bounds500, which(bounds500<Inf),500)
# create vector where all land above 100 yr flood is NA
bounds100 <- getValues(returnLevels_wse_upscaled[[which(return_periods==100)]])
# replace all non NA values with 100
bounds <- replace(bounds500, which(bounds100<Inf),100)
# set raster values to be either NA, 500, 100, or 1 to define the FEMA flood zones and mean WSE

# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas outside the 500 wet
vals<-replace(vals, which(is.na(getValues(floodzones))==TRUE), NA)
# set raster values to the interpolated return period values

if(map_type=='return period'){
  floodzones <- setValues(floodzones, bounds)
  rt_map <- setValues(rt_map, vals)
}
else if(map_type=='probability'){
  floodzones <- setValues(floodzones, 1/bounds)
  rt_map <- setValues(rt_map, 1/vals)
}
else if(map_type=='percent probability'){
  floodzones <- setValues(floodzones, 1/bounds*100)
  rt_map <- setValues(rt_map, 1/vals*100)
}
else(print("ERROR: define map type ('return period', 'probability', 'percent probability)"))


names(rt_map)<-'FLOPIT flood probability map'
names(floodzones)<-'FLOPIT flood zone map'
######################### Save analysis data, flood probability map, and flood zones map
if(save == TRUE){
# save the RDATA for analysis
save.image(file = save_path_data)

# save the rt_map (interpolated probability map using all FEMA flood surface elevation values)
writeRaster(rt_map, filename = save_path_map, format = 'GTiff', overwrite = TRUE)

# save the FEMA flood zone map (Warning: not guaranteed to exactly match the official FEMA NFHL flood zone map)
writeRaster(floodzones, filename = save_path_zones, format = 'GTiff', overwrite = TRUE)
}
######################### Return FLOPIT results (rt map and flood zone map)
if(save == FALSE){
  results <- list(rt_map, floodzones)
  return(results)
}
}
