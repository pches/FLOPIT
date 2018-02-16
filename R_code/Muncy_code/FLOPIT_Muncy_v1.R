###################################################
# file: Flood_Map_Interpolation_Clean_Muncy.R
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
# Last changes: February 16, 2018 (K. Joel Roop-Eckart)
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
#install.packages('scales')
require('scales')

# set working directory to a folder contain the R_code, Figures, and Data folders
setwd('....')

################ Import Houston Clipped and Resampled Raster files ###########################
# be sure to change file paths

# import 10% annual chance flood WSE raster
str_name <- 'Data/Muncy_data/depth_10pct_projected_clipped.tif'
depth_10pct <- raster(str_name)

# import 2% annual chance flood WSE raster
str_name <- 'Data/Muncy_data/depth_10pct_projected_clipped.tif'
depth_2pct <- raster(str_name)

# import 1% annual chance flood WSE raster
str_name <- 'Data/Muncy_data/depth_10pct_projected_clipped.tif'
depth_1pct <- raster(str_name)

# import 0.2% annual chance flood WSE raster
str_name <- 'Data/Muncy_data/depth_10pct_projected_clipped.tif'
depth_02pct <- raster(str_name)

# import Lidar DEM of the study area of Houston
str_name <- 'Data/Muncy_data/depth_10pct_projected_clipped.tif'
elev <- raster(str_name)


############ Convert flood depths (ft) to elevations (ft) ###############

# convert depths to elevations
vals_10pct <- getValues(depth_10pct)+getValues(elev)
wse_10pct <- setValues(depth_10pct, vals_10pct)

vals_2pct <- getValues(depth_2pct)+getValues(elev)
wse_2pct <- setValues(depth_10pct, vals_2pct)

vals_1pct <- getValues(depth_1pct)+getValues(elev)
wse_1pct <- setValues(depth_10pct, vals_1pct)

vals_02pct <- getValues(depth_02pct)+getValues(elev)
wse_02pct <- setValues(depth_10pct, vals_02pct)

################ Aggregate rasters to lower resolutions #####################
# rasters are too large for quick computation, aggregate by 20 times pixel width 
# (from ~16.4x16.4 to ~ 328x328 feet cells)
aggval <-20 # define aggregation size
# aggregate the rasters
wse_10pct_agg<-aggregate(wse_10pct, fact = aggval, fun = mean)
wse_2pct_agg<-aggregate(wse_2pct, fact = aggval, fun = mean)
wse_1pct_agg<-aggregate(wse_1pct, fact = aggval, fun = mean)
wse_02pct_agg<-aggregate(wse_02pct, fact = aggval, fun = mean)
elev_agg<-aggregate(elev, fact = aggval, fun = mean)

# calculate new raster dimensions after aggregation
rast_x <- elev_agg@ncols
rast_y <- elev_agg@nrows

# create matrix of points located at the center of each raster pixel for newly aggregated rasters
pixel <- (wse_02pct_agg@extent@xmax-elev_agg@extent@xmin)/elev_agg@ncols
xy <- matrix(data = NA, nrow = rast_x*rast_y,ncol = 2)
xy[,2] <- rep(seq((elev_agg@extent@xmin+0.4*pixel),(elev_agg@extent@xmax-0.4*pixel),pixel), each = rast_y)
xy[,1] <- rep(seq((elev_agg@extent@ymin+0.4*pixel),((elev_agg@extent@ymax-0.4*pixel)),pixel), rast_x)

######### Use Voronoi Polygon Interpolation/extrapolation ##########

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

# create voronoi polygons of wse_2pct for example figure
wse_2pct_points <- rasterToPoints(wse_2pct_agg, spatial = TRUE)
v2<-voronoi(wse_2pct_points)

# start spatial WSE interpolation/extrapolation
wse_10pct_interp <- poly_interp(wse_10pct_agg)
wse_2pct_interp <- poly_interp(wse_2pct_agg)
wse_1pct_interp <- poly_interp(wse_1pct_agg)
wse_02pct_interp <- poly_interp(wse_02pct_agg)

#################### Interpolate flood Return Period ###################
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
  floods <- sort(floods, decreasing = FALSE) # due to occasional error introduced by aggregation 
  # on steep slopes, smaller floods may produce higher WSE elevations for a single tile.
  # While these errors reduce accuracy, sorting the floods allows the interpolation scheme to make a 
  # realistic approximation or the relationship and finish the job
  vals[i] <- spline(floods,returns,xout = elev_agg@data@values[i], method = 'hyman')$y
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

# recreate FEMA flood boundaries
# create new raster
FEMA_floodbounds <- wse_02pct_agg
# create vector where all land above 500 yr flood is NA
bounds500 <- getValues(wse_02pct_agg)
# replace all non NA values with 500
bounds500 <- replace(bounds500, which(bounds500<10000),500)
# create vector where all land above 100 yr flood is NA
bounds100 <- getValues(wse_1pct_agg)
# replace all non NA values with 100
bounds100 <- replace(bounds500, which(bounds100<10000),100)
bounds <- bounds100
# set raster values to be either NA, 500, 100, or 1 to define the FEMA flood zones and mean WSE
FEMA_floodbounds <- setValues(FEMA_floodbounds, bounds100)

# set raster values to the interpolated return period values
rt_map <- setValues(rt_map, vals)

############################ Create Raster Map Figures ############################

# load colorblind approved colors for figures from colorbrewer
cols <- brewer.pal(10, 'Spectral')

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse10pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_10pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(480, 570, 10)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse2pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_2pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(480, 570, 10)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse1pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_1pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(480, 570, 10)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse02pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_02pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(480, 800, 10)))

dev.off()
pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_elev.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(elev_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(450, 650, 20)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_probmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_FEMA_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(FEMA_floodbounds, xlab = "Easting", ylab = "Northing",col.regions = cols,
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_voronoi_polygons.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(v2, xlab = "Easting", ylab = "Northing",
       par.settings = list(fontsize = list(text = 20))))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_wse10pct_forVoronoi.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_10pct_agg, xlab = "Easting", ylab = "Northing",
       par.settings = list(fontsize = list(text = 20))))
dev.off()

rt_map_error <- rt_map
rt_map_error <- setValues(rt_map_error, getValues(FEMA_floodbounds)-getValues(rt_map))

pdf("Figures/Muncy_figures/PCHES_flood_interpolation_probmap_error.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map_error, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 50)))
dev.off()

########### Create interpolation diagnostic figures #############

FEMA_rt <- bounds
interp_rt <- vals

# determine which cells posses NA values in any of the flood return period maps
allnas <- c(which(is.nan(FEMA_rt)==TRUE), which(is.nan(interp_rt)==TRUE),
            which(is.na(FEMA_rt)==TRUE), which(is.na(interp_rt)==TRUE))
length(allnas)
allnas <- unique.numeric_version(allnas)
length(allnas)

# place cell values with real values across all three return period maps into matrix
flood_dat <- matrix(data = NA, nrow = length(bounds[-allnas]), ncol = 2)
flood_dat[,1] <- bounds[-allnas]
flood_dat[,2] <- vals[-allnas]
colnames(flood_dat)<-as.character(c("FEMA_zones", "Interpolated_RTs"))

# plot FEMA flood zone return periods vs FEMA calculated and interpolated return periods 
pdf("Figures/Muncy_figures/PCHES_zone_vs_interp_boxplot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(Interpolated_RTs~FEMA_zones, data = flood_dat, col = 'lightblue',
        xlab = "FEMA Flood Zones (Return Period)", ylab = "Interpolated Return Period", pch = 16)
dev.off()

cols <- brewer.pal(3, 'Dark2')
pdf("Figures/Muncy_figures/PCHES_zone_vs_interp_scatterplot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
plot(flood_dat[,2], flood_dat[,1], ylab = "FEMA Flood Zone Return Period", 
     xlab = "Interpolated Return Period", xlim = c(0, 500), ylim = c(0, 500), pch = 16, col = alpha('black', 0.5))
points(1:500, 1:500, type = 'l', lwd = 3, col = 'blue')
points(1:100,rep(100, 100), type = 'l', lwd = 3, col = 'red')
points(101:500,rep(500, 400), type = 'l', lwd = 3, col = 'red')
points(rep(100, 400),101:500, type = 'l', lwd = 3, col = 'red')
legend('topleft', legend = c("Raster Cells", "1 to 1 Match", 'Flood Zones'), 
       col = c('black', 'blue', 'red'), pch = c(16, NA, NA), lwd = c(NA, 2, 2), bg = 'white')
dev.off()

############ Goodness of fit analysis and associated figures ########################
# 1: of all overlapping cells with flood probabilities for each scheme
# what's the average flood probability?

# FEMA flood map
print(c(round(mean(flood_dat[,1])), 'Average FEMA flood zone return period')) # FEMA flood zones: 123 yrs
print(c(round(mean(flood_dat[,2])), 'Average Interpolated flood return period')) # Interpolated flood return periods: 58 yrs

pdf("figures/PCHES_Flood_Interpolation_Histogram.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
hist(flood_dat[,1], col = cols[1], breaks = 20, density = 10, freq = FALSE, main = '')
hist(flood_dat[,2], col = cols[2], add = TRUE, breaks = 20, density = 10, freq = FALSE)
legend('topright', legend = c("FEMA Zones", "Interpolated Return Periods"), 
       col = c(cols[1],cols[2]), 
       pch = c(15, 15), lwd = c(NA, NA))
dev.off()

# 2: how many cells were errors of ommission vs commission for the probability maps vs FEMA flood zones

# Interpolated return period map errors of ommission
omitted_vals<- FEMA_rt[(intersect(which(is.na(interp_rt)==TRUE),which(is.na(FEMA_rt)==FALSE)))]
length(omitted_vals) # 196 errors of omission out of 2821
sum(omitted_vals==100) # 75 omissions in the 100yr zone
sum(omitted_vals==500) # 121 omissions in the 100yr zone
sum(omitted_vals==10) # 0 omissions in the 1yr zone (general WSE)

# Interpolated return period map errors of commission
committed_vals<- interp_rt[(intersect(which(is.na(interp_rt)==FALSE),which(is.na(FEMA_rt)==TRUE)))]
length(committed_vals) # 22 errors of commission
sum(committed_vals>100) # 22 commissions in the 100yr zone
sum(committed_vals<=100) # 0 commissions in the 100yr zone
sum(committed_vals==0) # 0 commissions in the 1yr zone (general WSE)

########################### Root Mean Square Error and associated figures ###################################

# additional error introduced by omission of values
arena1 <- c(unique.numeric_version(c(which(is.nan(interp_rt)==TRUE), which(is.na(interp_rt)==TRUE))))
notna1 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==FALSE & is.na(FEMA_rt)==FALSE))))
# list of values interpolation left out
# additional error introduced by commission of values
arena2 <- c(unique.numeric_version(c(which(is.nan(FEMA_rt)==TRUE), which(is.na(FEMA_rt)==TRUE))))
notna2 <- c(unique.numeric_version(c(which(is.nan(interp_rt)==FALSE & is.na(interp_rt)==FALSE))))
# list of values interpolation should not have included

# root mean square error, accounting for omission and commission errors
rmse_interp <- sqrt(mean((c(((flood_dat[,2]-flood_dat[,1])^2),
                            ((FEMA_rt[intersect(arena1,notna1)]-500)^2),
                            ((interp_rt[intersect(arena2,notna2)]-500)^2)))))
print(c(rmse_interp, 'RMSE'))

cdfvals <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfvals[i]<-((length(which(flood_dat[,2]<=i)))/(length(flood_dat[,2])))}

cdfbounds <- vector(mode = 'numeric', length = length(1:500))
for(i in 1:500){cdfbounds[i]<-((length(which(flood_dat[,1]<=i)))/(length(flood_dat[,1])))}

cols <- brewer.pal(3, 'Dark2')
pdf("Figures/Muncy_figures/PCHES_Muncy_generalhazard_plot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(mfrow=c(2,2))
# plot pdfs
hist(vals, freq = FALSE, col = cols[1], ylim = c(0,0.05), xlim = c(0,500), breaks = 20, density = 0, 
     main = '', xlab = 'Return period (yrs)', ylab = 'Density')
par(new = TRUE)
hist(bounds, freq = FALSE, col = cols[2], ylim = c(0,0.05), xlim = c(0,500), breaks = 20, density = 0,
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
par(new=TRUE)
plot(density(vals, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[1], lwd = 2, 
     main = '', xlab = '', ylab = '', xaxt='n',yaxt='n')
lines(density(bounds, na.rm = TRUE), ylim = c(0,0.05), xlim = c(0,500), col = cols[2], lwd = 2)
legend('topright', legend = c("Interpolated", "FEMA"), 
       col = c(cols[1],cols[2]), 
       pch = c(NA, NA), lwd = c(2,2))
# plot cdfs
plot(1:500,cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = 'Cumulative density')
lines(cdfbounds, col = cols[2], lwd = 2)
legend('bottomright', legend = c("Interpolated", "FEMA"), 
       col = c(cols[1],cols[2]), 
       pch = c(NA, NA), lwd = c(2,2))
# plot survival function
plot(1:500,1-cdfvals, type = 'l', col = cols[1], lwd = 2, xlab = 'Return period (yrs)',
     ylab = '1-cdf (cumulative density)')
lines(1-cdfbounds, col = cols[2], lwd = 2)
legend('topright', legend = c("Interpolated", "FEMA"), 
       col = c(cols[1],cols[2]), 
       pch = c(NA, NA), lwd = c(2,2))
plot(density(-c(((flood_dat[,2]-flood_dat[,1])),
               ((FEMA_rt[intersect(arena1,notna1)]-500)),
               ((interp_rt[intersect(arena2,notna2)]-500)))), xlab = 'Return period error (yrs)', 
     main = '', ylab = 'Error Density')
abline(v=0, col = 'green', lwd =2)
legend('topleft', legend= c('Error density', '0 error line'), lwd = c(1,2), col = c('black', 'green'), bg = 'white')
dev.off()
