###################################################
# file: Flood_Map_Interpolation_Clean_Houston.R
#
###################################################
# Author and copyright: K. Joel Roop-Eckart
# The Pennsylvania State University
# kjr30@psu.edu
#
# Distributed under the GNU general public license
# No warranty
#
###################################################
# Last changes: February 15, 2018 (K. Joel Roop-Eckart)
###################################################

rm(list = ls())
#dev.off()

# load necessary packages
#install.packages('raster')
require('raster')
#install.packages('rgdal')
require('rgdal')
#install.packages('sp')
require('sp')
#install.packages('FNN')
require('FNN')
#install.packages('dismo')
require('dismo')
#install.packages('deldir')
require('deldir')
#install.packages('rgeos')
require('rgeos')
#install.packages('RColorBrewer')
require('RColorBrewer')

# set working directory
setwd('C:/Users/kjr30/Desktop')
################ Import Houston Clipped and Resampled Raster files ###########################
# be sure to change file paths

# import 10% annual chance flood WSE raster
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/wse_10pct_clipped_resampled.tif'
wse_10pct <- raster(str_name)

# import 4% annual chance flood WSE raster 
# (NOTE: Due to data quality concerns, this raster is not used)

# import 2% annual chance flood WSE raster
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/wse_2pct_clipped_resampled.tif'
wse_2pct <- raster(str_name)

# import 1% annual chance flood WSE raster
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/wse_1pct_clipped_resampled.tif'
wse_1pct <- raster(str_name)

# import 0.2% annual chance flood WSE raster
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/wse_02pct_clipped_resampled.tif'
wse_02pct <- raster(str_name)

# import Lidar DEM of the study area of Houston
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/ned19_n29x75_w095x50_tx_houstoncity_2008_clipped_feet.tif'
elev <- raster(str_name)

# import FEMA calculated annual percent flood chance raster
str_name <- 'C:/Users/kjr30/Desktop/Data/Houston_data/pctannchance_clipped_resampled.tif'
prob_map_orig <- raster(str_name)

################ Extrapolating to missing values #####################

# rasters are too large for quick computation, aggregate by 5 times pixel width 
# (from ~10x10 to ~ 50x50 feet)
aggval <-5 # define aggregation size
# aggregate the rasters
wse_10pct_agg<-aggregate(wse_10pct, fact = aggval)
wse_2pct_agg<-aggregate(wse_2pct, fact = aggval)
wse_1pct_agg<-aggregate(wse_1pct, fact = aggval)
wse_02pct_agg<-aggregate(wse_02pct, fact = aggval)
elev_agg<-aggregate(elev, fact = aggval)
prob_map_orig_agg <- aggregate(prob_map_orig, fact = aggval)

# calculate new raster dimensions after aggregation
rast_x <- wse_10pct_agg@ncols
rast_y <- wse_10pct_agg@nrows

# create matrix of points located at the center of each raster pixel for newly aggregated rasters
pixel <- (wse_02pct_agg@extent@xmax-wse_02pct_agg@extent@xmin)/wse_02pct_agg@ncols
xy <- matrix(data = NA, nrow = rast_x*rast_y,ncol = 2)
xy[,2] <- rep(seq((wse_02pct_agg@extent@xmin+0.4*pixel),(wse_02pct_agg@extent@xmax-0.4*pixel),pixel), each = rast_y)
xy[,1] <- rep(seq((wse_02pct_agg@extent@ymin+0.4*pixel),((wse_02pct_agg@extent@ymax-0.4*pixel)),pixel), rast_x)

######### Use Voronoi Polygon Interpolation/extrapolation

# define voronoi interpolation/extrapolation function
poly_interp <- function(rast){
  rast_points <- rasterToPoints(rast, spatial = TRUE)
  v<-voronoi(rast_points)
  v_rast <- rasterize(v, rast, field = v@data)
  return(v_rast)
}

# create voronoi polygons of wse_10pct for example figure
wse_10pct_points <- rasterToPoints(wse_10pct_agg, spatial = TRUE)
v<-voronoi(wse_10pct_points)
# start spatial WSE interpolation/extrapolation
wse_10pct_interp <- poly_interp(wse_10pct_agg)
wse_2pct_interp <- poly_interp(wse_2pct_agg)
wse_1pct_interp <- poly_interp(wse_1pct_agg)
wse_02pct_interp <- poly_interp(wse_02pct_agg)

#################### Interpolate flood Return Period ###################
start<-Sys.time()
# calculate return period values for each flood
returns <- c(1/(10/100),1/(2/100),1/(1/100),1/(0.2/100)) # 4% flood and mean wse were removed due to data quality
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map <- wse_2pct_interp
# set all values of the new raster to be NA
setValues(rt_map, NA)

# define a vector vals to write interpolated flood return periods to
vals <- vector(mode = 'numeric', length = length(rt_map@data@values))
# for loop uses splines to define flood elevation to return period for each cell 
# and interpolate the return period associated with each cell's elevation
for (i in 1:length(rt_map@data@values)) {
  floods <- c(wse_10pct_interp@data@values[i],wse_2pct_interp@data@values[i],wse_1pct_interp@data@values[i],wse_02pct_interp@data@values[i])
    vals[i] <- spline(floods,returns,xout = elev_agg[i], method = 'hyman')$y
}

# return periods cannot be negative, replace negative return period with 0
vals<-replace(vals,which(vals<min(returns)),min(returns))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals<-replace(vals,which(vals>500),NA)

# "hyman" method only gurantees monotonically increasing interpolation for all values within the 
# data range. Replace all elevation values below the data range with lowest return period
for(i in 1:length(vals)){
  vals[i]<-replace(vals[i],which(min(wse_10pct_interp[i],wse_2pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])>elev_agg[i]),min(returns))
}
# replace all elevation values above the data range with NA
for(i in 1:length(vals)){
  vals[i]<-replace(vals[i],which(max(wse_10pct_interp[i],wse_2pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])<elev_agg[i]),NA)
}
end<-Sys.time()
splinespeed<-(end-start)
#################### Loess Case ###################
start<-Sys.time()

# calculate return period values for each flood
returns <- c(1/(10/100),1/(2/100),1/(1/100),1/(0.2/100)) # 4% flood and mean wse were removed due to data quality
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map_loess <- wse_2pct_interp
# set all values of the new raster to be NA
setValues(rt_map_loess, NA)

# define a vector vals to write interpolated flood return periods to
vals_loess <- vector(mode = 'numeric', length = length(rt_map_loess@data@values))
# for loop uses splines to define flood elevation to return period for each cell 
# and interpolate the return period associated with each cell's elevation
for (i in 1:length(rt_map_loess@data@values)) {
  floods <- c(wse_10pct_interp@data@values[i],wse_2pct_interp@data@values[i],wse_1pct_interp@data@values[i],wse_02pct_interp@data@values[i])
  fit <- loess(returns~floods, span = 1)
  vals_loess[i] <- predict(fit, elev_agg[i])
}

# return periods cannot be negative, replace negative return period with 0
vals_loess<-replace(vals_loess,which(vals_loess<min(returns)),min(returns))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals_loess<-replace(vals_loess,which(vals_loess>500),NA)

# "hyman" method only gurantees monotonically increasing interpolation for all values within the 
# data range. Replace all elevation values below the data range with lowest return period
for(i in 1:length(vals_loess)){
  vals_loess[i]<-replace(vals_loess[i],which(min(wse_10pct_interp[i],wse_2pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])>elev_agg[i]),min(returns))
}
# replace all elevation values above the data range with NA
for(i in 1:length(vals_loess)){
  vals_loess[i]<-replace(vals_loess[i],which(max(wse_10pct_interp[i],wse_2pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])<elev_agg[i]),NA)
}
end<-Sys.time()
loessspeed<-(end-start)
print(c(as.numeric(loessspeed)/as.numeric(splinespeed), "Spline interpolation performed ___ times faster than loess"))
############################ Missing data Case ########################################
# calculate return period values for each flood
returns_miss <- c(1/(10/100),1/(1/100),1/(0.2/100)) # 4% flood and mean wse were removed due to data quality
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map_miss <- wse_2pct_interp
# set all values of the new raster to be NA
setValues(rt_map_miss, NA)

# define a vector vals to write interpolated flood return periods to
vals_miss <- vector(mode = 'numeric', length = length(rt_map_miss@data@values))
# for loop uses splines to define flood elevation to return period for each cell 
# and interpolate the return period associated with each cell's elevation
for (i in 1:length(rt_map_miss@data@values)) {
  floods_miss <- c(wse_10pct_interp@data@values[i],wse_1pct_interp@data@values[i],wse_02pct_interp@data@values[i])
  vals_miss[i] <- spline(floods_miss,returns_miss,xout = elev_agg[i], method = 'hyman')$y
}

# return periods cannot be negative, replace negative return period with 0
vals_miss<-replace(vals_miss,which(vals_miss<min(returns_miss)),min(returns_miss))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals_miss<-replace(vals_miss,which(vals_miss>500),NA)

# "hyman" method only gurantees monotonically increasing interpolation for all values within the 
# data range. Replace all elevation values below the data range with lowest return period
for(i in 1:length(vals_miss)){
  vals_miss[i]<-replace(vals_miss[i],which(min(wse_10pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])>elev_agg[i]),min(returns_miss))
}
# replace all elevation values above the data range with NA
for(i in 1:length(vals_miss)){
  vals_miss[i]<-replace(vals_miss[i],which(max(wse_10pct_interp[i],wse_1pct_interp[i],wse_02pct_interp[i])<elev_agg[i]),NA)
}

##################### Create flood probability maps ################

# create a vector of cell flooding probabilities from the FEMA flood probability raster
probabilities<-getValues(prob_map_orig_agg)
# theoretically they could not calculate probabilities less than 0.002, so overwrite all probabilities
# less than 0.2% with NA.
probabilities<-replace(probabilities, which(probabilities<0.002),NA)
# convert probabilities to return periods
returnperiods <- 1/probabilities
# create new raster
rt_map_orig <- prob_map_orig_agg
# overwrite raster values with return periods
rt_map_orig <- setValues(rt_map_orig, returnperiods)

# recreate FEMA flood boundaries
# create new raster
FEMA_floodbounds <- wse_02pct_agg
# create vector where all land above 500 yr flood is NA
bounds500 <- getValues(wse_02pct_agg)
# replace all non NA values with 500
bounds500 <- replace(bounds500, which(bounds500<500),500)
# create vector where all land above 100 yr flood is NA
bounds100 <- getValues(wse_1pct_agg)
# replace all non NA values with 100
bounds <- replace(bounds500, which(bounds100<100),100)
# create vector of elevation values
boundsmean <- getValues(elev_agg)
# replace all values under 10 ft with 1 (visual estimate of assumed mean WSE)
bounds <- replace(bounds, which(boundsmean<1),1)
# set raster values to be either NA, 500, 100, or 1 to define the FEMA flood zones and mean WSE
FEMA_floodbounds <- setValues(FEMA_floodbounds, bounds)


# set raster values to the interpolated return period values
rt_map <- setValues(rt_map, vals)

# set raster values to the interpolated return period values
rt_map_loess <- setValues(rt_map_loess, vals_loess)

# set raster values to the interpolated return period values
rt_map_miss <- setValues(rt_map_miss, vals_miss)

############################ Plots For Presentation ############################
# be sure to change all file pathing for saved files

# load colorblind approved colors for figures from colorbrewer
cols <- brewer.pal(10, 'Spectral')

pdf("Figures/Houston_figures/PCHES_flood_interpolation_wse10pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_10pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(-1, 40, 4)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_wse2pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_2pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(-1, 40, 4)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_wse1pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_1pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(-1, 40, 4)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_wse02pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_02pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(-1, 40, 4)))

dev.off()
pdf("Figures/Houston_figures/PCHES_flood_interpolation_elev.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(elev_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(-1, 40, 4)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_probmap_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_miss_probmap_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_miss, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_loess_probmap_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_loess, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_probmap_orig.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_orig, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_FEMA_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(FEMA_floodbounds, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_voronoi_polygons.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(v, xlab = "Easting", ylab = "Northing",
       par.settings = list(fontsize = list(text = 20))))
dev.off()

pdf("Figures/Houston_figures/PCHES_wse10pct_forVoronoi.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_10pct_agg, xlab = "Easting", ylab = "Northing",
       par.settings = list(fontsize = list(text = 20))))
dev.off()


# create error map plots
  #calculate errors
rt_map_error <- rt_map
rt_map_error <- setValues(rt_map_error, getValues(rt_map)-getValues(rt_map_orig))
rt_map_miss_error <- rt_map_miss
rt_map_miss_error <- setValues(rt_map_miss_error, getValues(rt_map_miss)-getValues(rt_map_orig))
rt_map_loess_error <- rt_map_loess
rt_map_loess_error <- setValues(rt_map_loess_error, getValues(rt_map_loess)-getValues(rt_map_orig))
FEMA_floodbounds_error <- FEMA_floodbounds
FEMA_floodbounds_error <- setValues(FEMA_floodbounds_error, getValues(FEMA_floodbounds)-getValues(rt_map_orig))

pdf("Figures/Houston_figures/PCHES_flood_interpolation_probmap_error.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_miss_probmap_error.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_miss_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_interpolation_loess_probmap_error.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_loess_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
dev.off()

pdf("Figures/Houston_figures/PCHES_flood_FEMA_zones_error.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(FEMA_floodbounds_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(-200, 450, 75)))
dev.off()

#pdf("figures/PCHES_probmap_error_plots.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
#par(mfrow=c(2,2))
#print(spplot(rt_map_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
#             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
#print(spplot(rt_map_miss_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
#             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
#print(spplot(rt_map_loess_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
#             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
#print(spplot(FEMA_floodbounds_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
#             par.settings = list(fontsize = list(text = 20)), at = seq(-200, 450, 75)))
#dev.off()

########### plot of how interpolated flood probabilities land in FEMA flood zones #############
# this is a work in progress, but you may find it interesting

# Figure 1 and 2: 
# Plot of interpolated flood probabilities vs FEMA flood zone
# Plot of interpolated flood probabilities vs FEMA calculated flood probabilities

FEMA_rt <- bounds
FEMA_prob_rt <- returnperiods
interp_rt <- vals
interp_rt_miss <- vals_miss
interp_rt_loess <- vals_loess

# determine which cells posses NA values in any of the flood return period maps
allnas <- c(which(is.nan(FEMA_rt)==TRUE), which(is.nan(FEMA_prob_rt)==TRUE), which(is.nan(interp_rt)==TRUE),
            which(is.na(FEMA_rt)==TRUE), which(is.na(FEMA_prob_rt)==TRUE), which(is.na(interp_rt)==TRUE))
length(allnas)
allnas <- unique.numeric_version(allnas)
length(allnas)

# place cell values with real values across all three return period maps into matrix
flood_dat <- matrix(data = NA, nrow = length(bounds[-allnas]), ncol = 3)
flood_dat[,1] <- bounds[-allnas]
flood_dat[,2] <- returnperiods[-allnas]
flood_dat[,3] <- vals[-allnas]
colnames(flood_dat)<-as.character(c("FEMA_zones", "FEMA_RTs", "Interpolated_RTs"))

# plot FEMA flood zone return periods vs FEMA calculated and interpolated return periods 

#
pdf("Figures/Houston_figures/PCHES_zone_vs_interp_boxplot_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(Interpolated_RTs~FEMA_zones, data = flood_dat, col = 'lightblue',
        xlab = "FEMA Flood Zones (Return Period)", ylab = "Interpolated Return Period", pch = 16)
dev.off()
#
pdf("Figures/Houston_figures/PCHES_zone_vs_FEMAcalc_boxplot_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(FEMA_RTs~FEMA_zones, data = flood_dat, col = 'lightblue',
        xlab = "FEMA Flood Zones (Return Period)", ylab = "FEMA Calculated Return Period", pch = 16)
dev.off()


# replace figure
pdf("Figures/Houston_figures/PCHES_zone_vs_interp_scatterplot_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
plot(flood_dat[,1], flood_dat[,3], xlab = "FEMA Flood Zone Return Period", 
     ylab = "Interpolated Return Period", xlim = c(0, 500), ylim = c(0, 500), pch = 16)
points(1:500, 1:500, type = 'l', lwd = 3, col = 'green')
legend('topleft', legend = c("Raster Cells", "1 to 1 Match"), col = c('black', 'green'), 
       pch = c(16, NA), lwd = c(NA, 2))
dev.off()

# plot FEMA calculated flood return periods vs interpolated flood return periods
pdf("Figures/Houston_figures/PCHES_FEMAcalc_vs_interp_scatterplot_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
plot(flood_dat[,2], flood_dat[,3], xlab = "FEMA Calculated Return Period", 
     ylab = "Interpolated Return Period", xlim = c(0, 500), ylim = c(0, 500), pch = 16)
points(1:500, 1:500, type = 'l', lwd = 3, col = 'green')
legend('topleft', legend = c("Raster Cells", "1 to 1 Match"), col = c('black', 'green'), 
       pch = c(16, NA), lwd = c(NA, 2))
dev.off()
# if the interpolation method performed in 
# this script recreates FEMA estimated probabilities, plotted points should land on the green 1 to 1 line
# if there is no bias in the itnerpolation method performed here, points should roughly follow the 
# 1 to 1 line.
pdf("Figures/Houston_figures/PCHES_FEMAcalc_vs_interp_smoothscatterplot_control.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
smoothScatter(flood_dat[,2], flood_dat[,3], xlab = "FEMA Calculated Return Period", 
              ylab = "Interpolated Return Period", xlim = c(0, 500), ylim = c(0, 500), pch = 16,
              colramp = colorRampPalette(c('white', 'red')), nrpoints = 100)
points(1:500, 1:500, type = 'l', lwd = 3, col = 'black')
legend('topleft', legend = c("Outlier Raster Cells", "Raster Cell Density", "1 to 1 Match"), 
       col = c('black', 'red','black'), 
       pch = c(16,15, NA), lwd = c(NA, NA, 2))
dev.off()
############ Goodness of fit analysis and more plots ########################
# 1: of all overlapping cells with flood probabilities for each scheme
# what's the average flood probability?

print(c(mean(bounds,na.rm=TRUE), "FEMA flood zones average return period (yrs)"))         # FEMA flood zones: 379 yrs
print(c(mean(returnperiods,na.rm=TRUE), "FEMA probability map average return period (yrs)"))  # FEMA calculated return periods: 219 yrs
print(c(mean(vals,na.rm=TRUE), "Interpolated average return period (yrs)"))           # Interpolated flood return periods: 226 yrs
print(c(mean(rt_map_miss[],na.rm=TRUE), "Interpolated with missing data average return period (yrs)"))  # missing data case? 236 yrs
print(c(mean(rt_map_loess[],na.rm=TRUE), "Interpolated with loess average return period (yrs)"))  # missing data case? 307 yrs

# 2: how many cells were errors of ommission vs commission for the probability maps vs FEMA flood zones
# FEMA return period map errors of ommission
omitted_vals<- FEMA_rt[(intersect(which(is.na(FEMA_prob_rt)==TRUE),which(is.na(FEMA_rt)==FALSE)))]
length(omitted_vals) # 222 errors of omission
sum(omitted_vals==100) # 4 omissions in the 100yr zone
sum(omitted_vals==500) # 218 omissions in the 100yr zone
sum(omitted_vals==1) # 0 omissions in the 1yr zone (general WSE)

# FEMA return period map errors of commission
committed_vals<- FEMA_prob_rt[(intersect(which(is.na(FEMA_prob_rt)==FALSE),which(is.na(FEMA_rt)==TRUE)))]
length(committed_vals) # 0 errors of commission
sum(committed_vals==100) # 0 commissions in the 100yr zone
sum(committed_vals==500) # 0 commissions in the 100yr zone
sum(committed_vals==1) # 0 commissions in the 1yr zone (general WSE)

# Interpolated return period map errors of ommission
omitted_vals<- FEMA_rt[(intersect(which(is.na(interp_rt)==TRUE),which(is.na(FEMA_rt)==FALSE)))]
length(omitted_vals) # 224 errors of omission
sum(omitted_vals==100) # 2 omissions in the 100yr zone
sum(omitted_vals==500) # 222 omissions in the 100yr zone
sum(omitted_vals==1) # 0 omissions in the 1yr zone (general WSE)

# Interpolated return period map errors of commission
committed_vals<- interp_rt[(intersect(which(is.na(interp_rt)==FALSE),which(is.na(FEMA_rt)==TRUE)))]
length(committed_vals) # 141 errors of commission
sum(committed_vals>100) # 141 commissions in the 100yr zone
sum(committed_vals<=100) # 0 commissions in the 100yr zone
sum(committed_vals==0) # 0 commissions in the 1yr zone (general WSE)

# Interpolated return period map errors of ommission loess case
omitted_vals<- FEMA_rt[(intersect(which(is.na(interp_rt_loess)==TRUE),which(is.na(FEMA_rt)==FALSE)))]
length(omitted_vals) # 1606 errors of omission
sum(omitted_vals==100) # 639 omissions in the 500yr zone
sum(omitted_vals==500) # 826 omissions in the 100yr zone
sum(omitted_vals==1) # 141 omissions in the 1yr zone (general WSE)

# Interpolated return period map errors of commission loess case
committed_vals<- interp_rt_loess[(intersect(which(is.na(interp_rt_loess)==FALSE),which(is.na(FEMA_rt)==TRUE)))]
length(committed_vals) # 105 errors of commission
sum(committed_vals>100) # 105 commissions in the 500yr zone
sum(committed_vals<=100) # 0 commissions in the 100yr zone
sum(committed_vals==0) # 0 commissions in the 1yr zone (general WSE)

# by ommission and commission, FEMA flood return period map is more accurate

# commission and omission are determined by the match between the DEM used for interpolation and 
# the DEM used for original flood modeling and interpolation scheme (spline vs loess, etc)

########################### Root Mean Square Error and more plots ###################################

# FEMA FLOOD ZONES

# additional error introduced by omission of values
arena3 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==TRUE), which(is.na(FEMA_rt)==TRUE))))
notna3 <- c(unique.numeric_version(c(which(is.nan(FEMA_prob_rt)==FALSE & is.na(FEMA_prob_rt)==FALSE))))
# list of values interpolation left out
#((FEMA_prob_rt[intersect(arena3,notna3)]-500)^2)
# additional error introduced by commission of values
arena4 <- c(unique.numeric_version(c(which(is.nan(FEMA_prob_rt)==TRUE), which(is.na(FEMA_prob_rt)==TRUE))))
notna4 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==FALSE & is.na(FEMA_rt)==FALSE))))
# list of values interpolation should not have included
#((FEMA_rt[intersect(arena4,notna4)]-500)^2)


rmse_floodzones <- sqrt(mean((c(((flood_dat[,2]-flood_dat[,1])^2),
                                ((FEMA_prob_rt[intersect(arena3,notna3)]-500)^2),
                                ((FEMA_rt[intersect(arena4,notna4)]-500)^2)))))
print(c(rmse_floodzones, "RMSE: FEMA Flood Zones"))

# Interpolated Flood Map

# additional error introduced by omission of values
arena1 <- c(unique.numeric_version(c(which(is.nan(interp_rt)==TRUE), which(is.na(interp_rt)==TRUE))))
notna1 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==FALSE & is.na(FEMA_rt)==FALSE))))
# list of values interpolation left out
# additional error introduced by commission of values
arena2 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==TRUE), which(is.na(FEMA_rt)==TRUE))))
notna2 <- c(unique.numeric_version(c(which(is.nan(interp_rt)==FALSE & is.na(interp_rt)==FALSE))))
# list of values interpolation should not have included

# root mean square error, accounting for omission and commission errors
rmse_interp <- sqrt(mean((c(((flood_dat[,2]-flood_dat[,3])^2),
                            ((FEMA_rt[intersect(arena1,notna1)]-500)^2),
                            ((interp_rt[intersect(arena2,notna2)]-500)^2)))))
print(c(rmse_interp, 'RMSE: Interpolated Return Periods'))

# Interpolated missing 50 year flood raster

# additional error introduced by omission of values
arena5 <- c(unique.numeric_version(c(which(is.nan(interp_rt_miss)==TRUE), which(is.na(interp_rt_miss)==TRUE))))
notna5 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==FALSE & is.na(FEMA_rt)==FALSE))))
# list of values interpolation left out
# additional error introduced by commission of values
arena6 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==TRUE), which(is.na(FEMA_rt)==TRUE))))
notna6 <- c(unique.numeric_version(c(which(is.nan(interp_rt_miss)==FALSE & is.na(interp_rt_miss)==FALSE))))
# list of values interpolation should not have included

# root mean square error, accounting for omission and commission errors for interpolation with missing data
rmse_interp_miss <- sqrt(mean((c(((flood_dat[,2]-interp_rt_miss[-allnas])^2),
                            ((FEMA_rt[intersect(arena5,notna5)]-500)^2),
                            ((interp_rt_miss[intersect(arena6,notna6)]-500)^2)))))
print(c(rmse_interp_miss, 'RMSE: Interpolated Return Periods, missing data'))

# Interpolated using loess fit

# additional error introduced by omission of values
arena7 <- c(unique.numeric_version(c(which(is.nan(interp_rt_loess)==TRUE), which(is.na(interp_rt_loess)==TRUE))))
notna7 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==FALSE & is.na(FEMA_rt)==FALSE))))
# list of values interpolation left out
# additional error introduced by commission of values
arena8 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==TRUE), which(is.na(FEMA_rt)==TRUE))))
notna8 <- c(unique.numeric_version(c(which(is.nan(interp_rt_loess)==FALSE & is.na(interp_rt_loess)==FALSE))))
# list of values interpolation should not have included

# root mean square error, accounting for omission and commission errors for interpolation with missing data
rmse_interp_loess <- sqrt(mean((c(((returnperiods[-unique(c(which(is.na(interp_rt_loess)==TRUE), 
                                 which(is.na(returnperiods)==TRUE)))]-
                                   interp_rt_loess[-unique(c(which(is.na(interp_rt_loess)==TRUE),
                                   which(is.na(returnperiods)==TRUE)))])^2),
                                 ((FEMA_rt[intersect(arena7,notna7)]-500)^2),
                                 ((interp_rt_loess[intersect(arena8,notna8)]-500)^2)))))
print(c(rmse_interp_loess, 'RMSE: Interpolated Return Periods, missing data'))


# calculate the cumulative distributions
cdfvals <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfvals[i]<-((length(which(vals[-which(is.na(vals)==TRUE)]<=i)))/(length(vals[-which(is.na(vals)==TRUE)])))}

cdfvals_miss <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfvals_miss[i]<-((length(which(vals_miss[-which(is.na(vals_miss)==TRUE)]<=i)))/(length(vals_miss[-which(is.na(vals_miss)==TRUE)])))}

cdfvals_loess <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfvals_loess[i]<-((length(which(vals_loess[-which(is.na(vals_loess)==TRUE)]<=i)))/(length(vals_loess[-which(is.na(vals_loess)==TRUE)])))}

cdfbounds <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfbounds[i]<-((length(which(bounds[-which(is.na(bounds)==TRUE)]<=i)))/(length(bounds[-which(is.na(bounds)==TRUE)])))}

cdffemart <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdffemart[i]<-((length(which(returnperiods[-which(is.na(returnperiods)==TRUE)]<=i)))/(length(returnperiods[-which(is.na(returnperiods)==TRUE)])))}


cols <- brewer.pal(3, 'Dark2')
pdf("Figures/Houston_figures/PCHES_Houston_generalhazard_plot.pdf",width = 8.9/2.54, height = 24/2.54)
par(mfrow=c(4,1))
# plot pdfs
hist(vals, freq = FALSE, col = cols[1], ylim = c(0,0.05), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = 'Return period (yrs)', ylab = 'Density')
par(new = TRUE)
hist(bounds, freq = FALSE, col = cols[2], ylim = c(0,0.05), xlim = c(0,500), breaks = 20, density = 0,
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
par(new=TRUE)
hist(returnperiods, freq = FALSE, col = cols[3], ylim = c(0,0.05), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = '', ylab = '')
par(new=TRUE)
plot(density(vals, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[1], lwd = 2, 
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
lines(density(bounds, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[2], lwd = 2)
lines(density(returnperiods, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[3], lwd = 2)
legend('topleft', legend = c("Interpolated", "FEMA Zones", 'FEMA Return Periods'), 
       col = c(cols[1],cols[2],cols[3]), 
       pch = c(NA, NA,NA), lwd = c(2,2,2), bg = 'white', cex = 1)
# plot cdfs
plot(1:500,cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = 'Cumulative density')
lines(cdfbounds, col = cols[2], lwd = 2)
lines(cdffemart, col = cols[3], lwd = 2)
legend('topleft', legend = c("Interpolated", "FEMA Zones", 'FEMA Return Periods'), 
       col = c(cols[1],cols[2],cols[3]), 
       pch = c(NA, NA, NA), lwd = c(2,2,2), bg = 'white', cex = 1)
# plot survival function
plot(1:500,1-cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = '1-cdf (cumulative density)')
lines(1-cdfbounds, col = cols[2], lwd = 2)
lines(1-cdffemart, col = cols[3], lwd = 2)
legend('bottomleft', legend = c("Interpolated", "FEMA Zones", 'FEMA Return Periods'), 
       col = c(cols[1],cols[2], cols[3]), 
       pch = c(NA, NA, NA), lwd = c(2,2,2), bg = 'white', cex = 1)
# plot error
plot(density(-c(((flood_dat[,2]-flood_dat[,3])),
               ((FEMA_rt[intersect(arena1,notna1)]-500)),
               ((interp_rt[intersect(arena2,notna2)]-500)))), xlab = 'Return period error (yrs)', 
     main = '', ylab = 'Error Density', col = cols[1], lwd = 2)
lines(density(-c(((flood_dat[,2]-interp_rt_miss[-allnas])),
                 ((FEMA_rt[intersect(arena1,notna1)]-500)),
                 ((interp_rt_miss[intersect(arena2,notna2)]-500)))), col = cols[3], lwd = 2)
lines(density(-c(((flood_dat[,2]-flood_dat[,1])),
                ((FEMA_prob_rt[intersect(arena3,notna3)]-500)),
                ((FEMA_rt[intersect(arena4,notna4)]-500)))), col = cols[2], lwd = 2)

abline(v=0, col = 'black', lwd =1)
legend('topright', legend= c('Interpolation scheme', 'FEMA Zones', 'Interpolation missing data', 
                            'Ideal interpolation'), lwd = c(2,2,2,1), col = c(cols[1], cols[2], 
                             cols[3],'black', 'green'), bg = 'white', cex = 1)
dev.off()


cols <- brewer.pal(5, 'Dark2')
pdf("Figures/Houston_figures/PCHES_Houston_generalhazard_modelcomp_plot.pdf",width = 8.9/2.54, height = 24/2.54)
par(mfrow=c(4,1))

# plot pdfs
hist(vals, freq = FALSE, col = cols[1], ylim = c(0,0.015), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = 'Return period (yrs)', ylab = 'Density')
par(new = TRUE)
hist(vals_miss, freq = FALSE, col = cols[4], ylim = c(0,0.015), xlim = c(0,500), breaks = 20, density = 0,
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
par(new=TRUE)
hist(returnperiods, freq = FALSE, col = cols[3], ylim = c(0,0.015), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = '', ylab = '')
par(new=TRUE)
hist(vals_loess, freq = FALSE, col = cols[5], ylim = c(0,0.015), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = '', ylab = '')
par(new=TRUE)
plot(density(vals, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[1], lwd = 2, 
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
lines(density(vals_miss, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[4], lwd = 2)
lines(density(vals_loess, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[5], lwd = 2)
lines(density(returnperiods, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[3], lwd = 2)
legend('topright', legend = c("Interpolated", "Interpolated (missing data)", "Interpolated (loess)", 
                              'FEMA Return Periods'), bg = 'white',
       col = c(cols[1],cols[5],cols[4],cols[3]), 
       pch = c(NA, NA,NA,NA,NA), lwd = c(2,2,2,2,2), cex = 1)

# plot cdfs
plot(1:500,cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = 'Cumulative density')
lines(cdfbounds, col = cols[2], lwd = 2)
lines(cdffemart, col = cols[3], lwd = 2)
legend('topright', legend = c("Interpolated", "Interpolated (missing data)", "Interpolated (loess)", 
                              'FEMA Return Periods'), bg = 'white',
       col = c(cols[1],cols[5],cols[4],cols[3]), 
       pch = c(NA, NA,NA,NA,NA), lwd = c(2,2,2,2,2), cex = 1)
# plot survival function
plot(1:500,1-cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = '1-cdf (cumulative density)')
#lines(1-cdfbounds, col = cols[2], lwd = 2)
lines(1-cdffemart, col = cols[3], lwd = 2)
lines(1-cdfvals_miss, col = cols[5], lwd = 2)
lines(1-cdfvals_loess, col = cols[4], lwd = 2)
legend('topright', legend = c("Interpolated", "Interpolated (missing data)", "Interpolated (loess)", 
                              'FEMA Return Periods'), bg = 'white',
       col = c(cols[1],cols[5],cols[4],cols[3]), 
       pch = c(NA, NA,NA,NA,NA), lwd = c(2,2,2,2,2), cex = 1)
# plot error
plot(density(-c(((flood_dat[,2]-flood_dat[,3])),
                ((FEMA_rt[intersect(arena1,notna1)]-500)),
                ((interp_rt[intersect(arena2,notna2)]-500)))), xlab = 'Return period error (yrs)', 
     main = '', ylab = 'Error Density', col = cols[1], lwd = 2)
lines(density(-c(((flood_dat[,2]-interp_rt_miss[-allnas])),
                 ((FEMA_rt[intersect(arena1,notna1)]-500)),
                 ((interp_rt_miss[intersect(arena2,notna2)]-500)))), col = cols[3], lwd = 2)
lines(density(-c(((returnperiods[-unique(c(which(is.na(interp_rt_loess)==TRUE), 
                                           which(is.na(returnperiods)==TRUE)))]-
                     interp_rt_loess[-unique(c(which(is.na(interp_rt_loess)==TRUE),
                                               which(is.na(returnperiods)==TRUE)))])),
                 ((FEMA_rt[intersect(arena7,notna7)]-500)),
                 ((interp_rt_loess[intersect(arena8,notna8)]-500)))), col = cols[5], lwd = 2)

abline(v=0, col = 'black', lwd =2)
legend('topright', legend = c("Interpolated", "Interpolated (missing data)", "Interpolated (loess)", 
                              'FEMA Return Periods'), bg = 'white',
       col = c(cols[1],cols[5],cols[4],cols[3]), 
       pch = c(NA, NA,NA,NA,NA), lwd = c(2,2,2,2,2), cex = 1)
dev.off()
