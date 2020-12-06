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
FLOPIT <- function(flood_rasters_names,flood_probabilities,elevation_raster,depth,aggregation_value,method,map_type,
                   save_outputs,save_path_data,save_path_map,save_path_zones){
  # Descriptions:
  # flood_rasters_names, 
  # flood_probabilities, 
  # elevation_raster, 
  # depth, 
  # aggregation_value, 
  # method, 
  # map_type,
  # save_outputs: Logical: If output data should be saved 
  # save_path_data: save_outputs 
  # save_path_map, 
  # save_path_zones
  
## For debugging the function, uncomment the following lines 
#flood_rasters_names=flood_rasters 
#flood_probabilities=flood_probabilities 
#elevation_raster=elevation_raster 
#depth = FALSE 
#aggregation_value = 1 
#method = 'log-linear'
#map_type = 'return period'
#save_outputs = FALSE 
  
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
dem <- raster(elevation_raster,values=TRUE)
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
if(aggregation_value > 1){
  returnLevels_wse_upscaled <- vector(length = n_rp)

  for(i in 1:n_rp){
    returnLevels_wse_upscaled[i] <- list(aggregate(flood_wse_vector[[i]], fact = aggregation_value, fun = mean))
  }

  dem_upscaled <- aggregate(dem, 
                      fact = aggregation_value, 
                      fun = mean)
}else{
  returnLevels_wse_upscaled <- flood_wse_vector
  dem_upscaled <- dem
}

print('Done with aggregating elevation and water surface elevation data')

######################### Smooth water surface elevation data using Inverse Weighted Distance method
# define inverse weighted distance interpolation/extrapolation function
# define raster grid that covers entire interpolation/extrapolation space (DEM suggested for this use)
flood_raster_grid <- rasterToPoints(dem_upscaled, spatial = TRUE)
gridded(flood_raster_grid) = TRUE

# define function
IWD_interp <- function(raster,grid){
    raster_points <- rasterToPoints(raster, spatial = TRUE)
    raster.gstat <- gstat(formula = rasterdata ~ 1, 
                          data = raster_points, 
                          nmax = 10, # for local kriging: the number of nearest observations that should be used for a kriging prediction or simulation
                          set = list(idp = 0.5)
                          )
    z <- predict(raster.gstat, grid)
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
flopit_interpolated_raster <- dem_upscaled

# set all values of the new raster to be NA
setValues(flopit_interpolated_raster, NA)

# create a vector to write interpolated flood return periods to (values are not assigned yet)
flopit_interpolated_rp_vec <- vector(mode = 'numeric', length = length(flopit_interpolated_raster@data@values))
# for loop uses splines to define flood elevation to return period for each cell
# and interpolate the return period associated with each cell's elevation
ground_elevations <- getValues(dem_upscaled)
n_points <- length(ground_elevations)


# For monitoring the progress
pb <- txtProgressBar(min = 0, max = length(getValues(flopit_interpolated_raster)), initial = 0, char = '=', style = 1)
start <- Sys.time()

for (i in 1:n_points){
  # monitor the progress 
  setTxtProgressBar(pb, i)
  # start the vector that contains the water surface evelations in this point 
  cell_wse <- vector(length = n_rp)
  for (j in 1:n_rp) {
    tryCatch(cell_wse[j] <- returnLevels_wse_smooth[[j]]@data@values[i]
             ,error=function() cell_wse[j] <- NA)
  }
  
  cell_wse <- sort(cell_wse, decreasing = FALSE) # due to occasional error introduced by aggregation, 
  # on steep slopes, smaller floods may produce higher WSE elevations for a single tile due to aggregation.
  # While these errors reduce accuracy, sorting the floods allows the interpolation scheme to make a
  # realistic approximation or the relationship and finish the job
  if(method == 'log-linear'){ # official FEMA log-linear interpolation formula
        denom1 <- 10^((
          (ground_elevations[i] - max(cell_wse[which(cell_wse < ground_elevations[i])]))*
            (log10(max(probs_pct[which(cell_wse > ground_elevations[i])])) - log10(min(probs_pct[which(cell_wse < ground_elevations[i])])))/
            (min(cell_wse[which(cell_wse > ground_elevations[i])]) - max(cell_wse[which(cell_wse < ground_elevations[i])]))
        ) + log10(min(probs_pct[which(cell_wse < ground_elevations[i])])))
        denom <- denom1/100 
        flopit_interpolated_rp_vec[i] <- 1/ denom
        
  }else if(method == 'spline'){
    
        flopit_interpolated_rp_vec[i] <- spline(cell_wse,return_periods,xout = ground_elevations[i], method = 'hyman')$y
  }

}

# finish progress monitoring
close(pb)
end <- Sys.time()
end-start

######################### Post process interpolated flood probability data
# return periods cannot be negative, replace negative return period with the minimum return period in the range
flopit_interpolated_rp_vec <- replace(
                                      flopit_interpolated_rp_vec,
                                      which(flopit_interpolated_rp_vec < min(return_periods)),
                                      min(return_periods))

# return periods over 500 are extrapolating beyond the data, set them to NA
flopit_interpolated_rp_vec <- replace(
                                      flopit_interpolated_rp_vec,
                                      which(flopit_interpolated_rp_vec > max(return_periods)),
                                      NA)

# if the wse of the largest flood is less than the elevation, set the pixel value to NA.
flopit_interpolated_rp_vec <- replace(
                                      flopit_interpolated_rp_vec, 
                                      which(getValues(returnLevels_wse_smooth[[n_rp]]) < getValues(dem_upscaled)), 
                                      NA)

# if the wse of the a flood is greater than the elevation, and the return period is NA, 
# set the pixel value to the return period of that flood.
for (i in 1:n_rp) {
      flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, 
                                        intersect(
                                                  which(getValues(returnLevels_wse_smooth[[i]]) > getValues(dem_upscaled)), 
                                                  which(is.na(flopit_interpolated_rp_vec) == TRUE)), 
                                        return_periods[i])
}

for (i in 1:n_rp) {
      flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, 
                                            intersect( 
                                                      which(is.na(getValues(returnLevels_wse_upscaled[[i]])) == FALSE), 
                                                      which(flopit_interpolated_rp_vec > return_periods[i])),
                                            return_periods[i])
}

flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, 
                                      intersect(
                                                which(is.na(getValues(returnLevels_wse_upscaled[[n_rp]])) == FALSE), 
                                                which(is.na(flopit_interpolated_rp_vec) == TRUE)), 
                                      return_periods[n_rp])

flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, 
                                      which(is.na(getValues(returnLevels_wse_upscaled[[n_rp]])) == TRUE), 
                                      NA)

# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in 
# the 500 year zone in the 100 year zone

for (i in 1:(n_rp-1)) {
  flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, intersect(intersect(which(is.na(getValues(returnLevels_wse_upscaled[[i+1]]))==FALSE), 
                                          which(is.na(getValues(returnLevels_wse_upscaled[[i]]))==TRUE)),
                                which((flopit_interpolated_rp_vec<return_periods[i]))), return_periods[i])
}

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
flopit_interpolated_rp_vec <- replace(flopit_interpolated_rp_vec, 
                                      which(is.na(getValues(floodzones))==TRUE), 
                                      NA)

# set raster values to the interpolated return period values
if(map_type == 'return period'){
  
      floodzones <- setValues(floodzones, bounds)
      flopit_interpolated_raster <- setValues(flopit_interpolated_raster, flopit_interpolated_rp_vec)
      
}else if(map_type == 'probability'){
  
      floodzones <- setValues(floodzones, 1/bounds)
      flopit_interpolated_raster <- setValues(flopit_interpolated_raster, 1/flopit_interpolated_rp_vec)
      
}else if(map_type=='percent probability'){
  
      floodzones <- setValues(floodzones, 1/bounds*100)
      flopit_interpolated_raster <- setValues(flopit_interpolated_raster, 1/flopit_interpolated_rp_vec*100)
      
}else{
  
      print("ERROR: define map type ('return period', 'probability', 'percent probability)")
}


names(flopit_interpolated_raster) <- 'FLOPIT flood probability map'
names(floodzones) <- 'FLOPIT flood zone map'

######################### Save analysis data, flood probability map, and flood zones map
if(save_outputs == TRUE){
  
      # save the RDATA for analysis
      save.image(file = save_path_data)

      # save the flopit_interpolated_raster (interpolated probability map using all FEMA flood surface elevation values)
      writeRaster(flopit_interpolated_raster, filename = save_path_map, format = 'GTiff', overwrite = TRUE)

      # save the FEMA flood zone map (Warning: not guaranteed to exactly match the official FEMA NFHL flood zone map)
      writeRaster(floodzones, filename = save_path_zones, format = 'GTiff', overwrite = TRUE)
      
}

results <- list(flopit_interpolated_raster, floodzones)
return(results)

} # End of FLOPIT function
