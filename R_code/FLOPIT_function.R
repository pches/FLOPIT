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
FLOPIT <- function(flood_rasters, flood_probabilities, elevation_raster, depth, aggregation_value, method, map_type,
                    save, save_path_data, save_path_map, save_path_zones){
################ Import Houston Clipped and Resampled Raster files ###########################
# be sure to change file paths
len <- length(flood_rasters)

flood_raster <- vector(length = len)
for(i in 1:len){
str_name <- flood_rasters[i]
flood_raster[i] <- list(raster(str_name))
names(flood_raster[[i]])<-'rasterdata'
}

# import Lidar DEM of the study area
str_name <- elevation_raster
elev <- raster(str_name)


######################### Convert flood depths (ft) to elevations (ft)
if(depth==TRUE)(
# convert depths to elevations
  for (i in 1:len) {
    rast_vals <- getValues(flood_raster[[i]])+getValues(elev)
    flood_raster[i] <- list(setValues(flood_raster[[i]], rast_vals))
  }
)
######################### Interpolate missing values
aggval <- aggregation_value # define aggregation size
# aggregate the rasters
flood_raster_agg <- vector(length = len)
for(i in 1:len){
  flood_raster_agg[i] <- list(aggregate(flood_raster[[i]], fact = aggval, fun = mean))
}

elev_agg<-aggregate(elev, fact = aggval, fun = mean)

######################### Use Inverse Weighted Distance Interpolation/extrapolation
# define inverse weighted distance interpolation/extrapolation function
# define raster grid that covers entire interpolation/extrapolation space (DEM suggested for this use)
rast_grid <- rasterToPoints(elev_agg, spatial = TRUE)
rast.grid <- rast_grid
gridded(rast.grid)=TRUE
# define function
IWD_interp <- function(rast,rast.grid){
  rast_points <- rasterToPoints(rast, spatial = TRUE)
  rast.gstat <- gstat(formula = rasterdata ~ 1, data = rast_points, 
                      nmax = 10, set = list(idp = 0.5))
  z <- predict(rast.gstat, rast.grid)
  #  z_rast <- rasterize(z, rast, field = z@data[,1])
  rast@data@values<-z@data[,1]
  return(rast)
}

# start spatial WSE interpolation/extrapolation
flood_raster_interp <- vector(length = len)
for (i in 1:len) {
  flood_raster_interp[i] <- list(IWD_interp(flood_raster_agg[[i]], rast.grid))
}

# Test for NA values after interpolation: If there are any, interpolation has made a mistake
for (i in 1:len) {
  out <- which(is.na(flood_raster_interp[[i]]@data@values)==TRUE)
  if(length(out)>0)(print('ERROR: Water surface elevation interpolation failed'))
}

######################### Interpolate flood probability map
# calculate return period values for each flood
returns <- 1/flood_probabilities # 4% flood and mean wse were removed due to data quality
probs <- 1/returns*100 # probabilities associated with the return periods
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map <- elev_agg
# set all values of the new raster to be NA
setValues(rt_map, NA)
# define a vector vals to write interpolated flood return periods to
vals <- vector(mode = 'numeric', length = length(rt_map@data@values))
# for loop uses splines to define flood elevation to return period for each cell
# and interpolate the return period associated with each cell's elevation
elevations <- getValues(elev_agg)
pb <- txtProgressBar(min = 0, max = length(getValues(rt_map)), initial = 0, char = '=', style = 1)
start <- Sys.time()
floods <- vector(length = len)
for (i in 1:length(rt_map@data@values)) {
  setTxtProgressBar(pb, i)
  for (j in 1:len) {floods[j] <- flood_raster_interp[[j]]@data@values[i]}
  floods <- sort(floods, decreasing = FALSE) # due to occasional error introduced by aggregation, 
  # on steep slopes, smaller floods may produce higher WSE elevations for a single tile due to aggregation.
  # While these errors reduce accuracy, sorting the floods allows the interpolation scheme to make a
  # realistic approximation or the relationship and finish the job
  if(method == 'log-linear')(vals[i] <- 1/((10^((
    (elevations[i]-max(floods[which(floods<elevations[i])]))*
      
      (log10(max(probs[which(floods>elevations[i])]))-log10(min(probs[which(floods<elevations[i])])))/
      
      (min(floods[which(floods>elevations[i])])-max(floods[which(floods<elevations[i])]))
    
  )+log10(min(probs[which(floods<elevations[i])]))))/100)) # official FEMA log-linear interpolation formula
  else if(method == 'spline')(vals[i] <- spline(floods,returns,xout = elevations[i], method = 'hyman')$y)
}
close(pb)
end <- Sys.time()
end-start

######################### Post process interpolated flood probability data
# return periods cannot be negative, replace negative return period with the minimum return period in the range
vals<-replace(vals,which(vals<min(returns)),min(returns))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals<-replace(vals,which(vals>max(returns)),NA)
# if the wse of the largest flood is less than the elevation, set the pixel value to NA.
vals<-replace(vals, which(getValues(flood_raster_interp[[len]])<getValues(elev_agg)), NA)
# if the wse of the a flood is greater than the elevation, and the return period is NA, 
# set the pixel value to the return period of that flood.

for (i in 1:len) {
  vals<-replace(vals, intersect(which(getValues(flood_raster_interp[[i]])>getValues(elev_agg)), which(is.na(vals)==TRUE)), returns[i])
}

for (i in 1:len) {
  vals<-replace(vals, intersect(which(is.na(getValues(flood_raster_agg[[i]]))==FALSE), which(vals>returns[i])), returns[i])
}

vals<-replace(vals, intersect(which(is.na(getValues(flood_raster_agg[[len]]))==FALSE), which(is.na(vals)==TRUE)), returns[len])

vals<-replace(vals, which(is.na(getValues(flood_raster_agg[[len]]))==TRUE), NA)
# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in 
# the 500 year zone in the 100 year zone

for (i in 1:(len-1)) {
  vals<-replace(vals, intersect(intersect(which(is.na(getValues(flood_raster_agg[[i+1]]))==FALSE), 
                                          which(is.na(getValues(flood_raster_agg[[i]]))==TRUE)),
                                which((vals<returns[i]))), returns[i])
}


#which(is.na(getValues(flood_raster_agg[[i]]))==TRUE)

########################## recreate flood zones
# create new raster
floodzones <- elev_agg
# create vector where all land above 500 yr flood is NA
bounds500 <- getValues(flood_raster_agg[[which(returns==500)]])
# replace all non NA values with 500
bounds500 <- replace(bounds500, which(bounds500<Inf),500)
# create vector where all land above 100 yr flood is NA
bounds100 <- getValues(flood_raster_agg[[which(returns==100)]])
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
