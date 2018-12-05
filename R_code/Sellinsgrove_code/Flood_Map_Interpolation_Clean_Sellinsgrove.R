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
# Last changes: February 14, 2018 (K. Joel Roop-Eckart)
###################################################

rm(list = ls())
#dev.off()

# Install packages
#install.packages(c('raster','rgdal','sp','FNN','dismo','deldir','rgeos','RColorBrewer', 'scales', 'gstat'))

# load necessary packages
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

########################## Links to project data

# Lidar
# http://maps.psiee.psu.edu/ImageryNavigator/

# FEMA data
# https://msc.fema.gov/portal/advanceSearch

# set working directory
setwd("~/GIS_Data/FLOPIT-master")
################ Import Houston Clipped and Resampled Raster files ###########################
# be sure to change file paths

# import 10% annual chance flood WSE raster
str_name <- 'Data/Sellinsgrove_data/dg_10pct_clip.tif'
depth_10pct <- raster(str_name)
names(depth_10pct)<-'rasterdata'

# import 4% annual chance flood WSE raster
str_name <- 'Data/Sellinsgrove_data/dg_4pct_clip.tif'
depth_4pct <- raster(str_name)
names(depth_10pct)<-'rasterdata'

# import 2% annual chance flood WSE raster
str_name <- 'Data/Sellinsgrove_data/dg_2pct_clip.tif'
depth_2pct <- raster(str_name)
names(depth_10pct)<-'rasterdata'

# import 1% annual chance flood WSE raster
str_name <- 'Data/Sellinsgrove_data/dg_1pct_clip.tif'
depth_1pct <- raster(str_name)
names(depth_10pct)<-'rasterdata'

# import 0.2% annual chance flood WSE raster
str_name <- 'Data/Sellinsgrove_data/dg_02pct_clip.tif'
depth_02pct <- raster(str_name)
names(depth_10pct)<-'rasterdata'

# import Lidar DEM of the study area of Houston
str_name <- 'Data/Sellinsgrove_data/Sellinsgrove_DEM_clip1.tif'
elev <- raster(str_name)
names(depth_10pct)<-'rasterdata'


############ Convert flood depths (ft) to elevations (ft) ###############

# convert depths to elevations
vals_10pct <- getValues(depth_10pct)+getValues(elev)
wse_10pct <- setValues(depth_10pct, vals_10pct)

vals_4pct <- getValues(depth_4pct)+getValues(elev)
wse_4pct <- setValues(depth_10pct, vals_4pct)

vals_2pct <- getValues(depth_2pct)+getValues(elev)
wse_2pct <- setValues(depth_10pct, vals_2pct)

vals_1pct <- getValues(depth_1pct)+getValues(elev)
wse_1pct <- setValues(depth_10pct, vals_1pct)

vals_02pct <- getValues(depth_02pct)+getValues(elev)
wse_02pct <- setValues(depth_10pct, vals_02pct)

################ Interpolate missing values #####################
# rasters are too large for quick computation, aggregate by 20 times pixel width
# (from ~16.4x16.4 to ~ 328x328 feet cells)
aggval <-10 # define aggregation size
# aggregate the rasters
wse_10pct_agg<-aggregate(wse_10pct, fact = aggval, fun = mean)
wse_4pct_agg<-aggregate(wse_4pct, fact = aggval, fun = mean)
wse_2pct_agg<-aggregate(wse_2pct, fact = aggval, fun = mean)
wse_1pct_agg<-aggregate(wse_1pct, fact = aggval, fun = mean)
wse_02pct_agg<-aggregate(wse_02pct, fact = aggval, fun = mean)
elev_agg<-aggregate(elev, fact = aggval, fun = mean)

################# Use Inverse Weighted Distance Interpolation/extrapolation ###################
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
start <- Sys.time()
wse_10pct_interp <- IWD_interp(wse_10pct_agg, rast.grid)
wse_4pct_interp <- IWD_interp(wse_4pct_agg, rast.grid)
wse_2pct_interp <- IWD_interp(wse_2pct_agg, rast.grid)
wse_1pct_interp <- IWD_interp(wse_1pct_agg, rast.grid)
wse_02pct_interp <- IWD_interp(wse_02pct_agg, rast.grid)
end <- Sys.time()
end-start
# time: ~ 16 second for five 2.8 million thousand pixel rasters

# Test for NA values after interpolation: If there are any, interpolation has gone wrong
which(is.na(wse_10pct_interp@data@values)==TRUE)
which(is.na(wse_2pct_interp@data@values)==TRUE)
which(is.na(wse_1pct_interp@data@values)==TRUE)
which(is.na(wse_02pct_interp@data@values)==TRUE)

#################### Interpolate flood Return Period ###################
# calculate return period values for each flood
returns <- c(1/(10/100),1/(4/100),1/(2/100),1/(1/100),1/(0.2/100)) # 4% flood and mean wse were removed due to data quality
# use pre-existing raster to create new raster of same cell size, extent, projection, etc.
rt_map <- wse_2pct_interp
# set all values of the new raster to be NA
setValues(rt_map, NA)

# define a vector vals to write interpolated flood return periods to
vals <- vector(mode = 'numeric', length = length(rt_map@data@values))
# for loop uses splines to define flood elevation to return period for each cell
# and interpolate the return period associated with each cell's elevation
elevations <- getValues(elev_agg)
pb <- txtProgressBar(min = 0, max = length(getValues(rt_map)), initial = 0, char = '=', style = 1)
for (i in 1:length(rt_map@data@values)) {
  setTxtProgressBar(pb, i)
  floods <- c(wse_10pct_interp@data@values[i],wse_4pct_interp@data@values[i],wse_2pct_interp@data@values[i],wse_1pct_interp@data@values[i],wse_02pct_interp@data@values[i])
  floods <- sort(floods, decreasing = FALSE) # due to occasional error introduced by aggregation
  # on steep slopes, smaller floods may produce higher WSE elevations for a single tile.
  # While these errors reduce accuracy, sorting the floods allows the interpolation scheme to make a
  # realistic approximation or the relationship and finish the job
  vals[i] <- spline(floods,returns,xout = elevations[i], method = 'hyman')$y
}
close(pb)
# return periods cannot be negative, replace negative return period with the minimum return period in the range
vals<-replace(vals,which(vals<min(returns)),min(returns))
# return periods over 500 are extrapolating beyond the data, set them to NA
vals<-replace(vals,which(vals>500),NA)
# if the wse of the largest flood is less than the elevation, set the pixel value to NA.
vals<-replace(vals, which(getValues(wse_02pct_interp)<getValues(elev_agg)), NA)
# if the wse of the a flood is greater than the elevation, and the return period is NA, 
# set the pixel value to the return period of that flood.
vals<-replace(vals, intersect(which(getValues(wse_10pct_interp)>getValues(elev_agg)), which(is.na(vals)==TRUE)), 10)
vals<-replace(vals, intersect(which(getValues(wse_2pct_interp)>getValues(elev_agg)), which(is.na(vals)==TRUE)), 50)
vals<-replace(vals, intersect(which(getValues(wse_1pct_interp)>getValues(elev_agg)), which(is.na(vals)==TRUE)), 100)
vals<-replace(vals, intersect(which(getValues(wse_02pct_interp)>getValues(elev_agg)), which(is.na(vals)==TRUE)), 500)

# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in the 100 year dry
vals<-replace(vals, intersect(which(is.na(getValues(wse_1pct_agg))==FALSE), which(vals>100)), 100)
# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in the 500 year dry
vals<-replace(vals, intersect(which(is.na(getValues(wse_02pct_agg))==FALSE), which(is.na(vals)==TRUE)), 500)

# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas in 
# the 500 year zone in the 100 year zone
vals<-replace(vals, intersect(intersect(which(is.na(getValues(wse_02pct_agg))==FALSE), 
                              which(is.na(getValues(wse_1pct_agg))==TRUE)),
              which((vals<100))), 100.1)

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
FEMA_floodbounds <- setValues(FEMA_floodbounds, bounds)

FEMA_floodbounds_prob <- setValues(FEMA_floodbounds, 1/bounds*100)


# Ensure proper flooding extent: Slight mistmatches between flood elevation and the DEM may leave areas outside the 500 wet
vals<-replace(vals, which(is.na(getValues(FEMA_floodbounds))==TRUE), NA)
# set raster values to the interpolated return period values
rt_map <- setValues(rt_map, vals)

prob_map <- setValues(rt_map, 1/vals*100)


############################ Save rt_map and FEMA flood bounds ##########################
# save the RDATA for analysis
save.image(file = "Data/Sellinsgrove_data/Sellinsgrove.RData")

# save the rt_map (interpolated probability map using all FEMA flood surface elevation values)
writeRaster(rt_map, filename = "Data/Sellinsgrove_data/FLOPIT_rt_map.tif", format = 'GTiff', overwrite = TRUE)

# save the rt_map (interpolated probability map using all FEMA flood surface elevation values)
writeRaster(FEMA_floodbounds, filename = "Data/Sellinsgrove_data/FEMA_zone_map.tif", format = 'GTiff', overwrite = TRUE)

# save the rt_map (interpolated probability map using all FEMA flood surface elevation values)
writeRaster(prob_map, filename = "Data/Sellinsgrove_data/FLOPIT_prob_map.tif", format = 'GTiff', overwrite = TRUE)

# save the rt_map (interpolated probability map using all FEMA flood surface elevation values)
writeRaster(FEMA_floodbounds_prob, filename = "Data/Sellinsgrove_data/FEMA_zone_prob_map.tif", format = 'GTiff', overwrite = TRUE)

############################ Plots ############################

# load colorblind approved colors for figures from colorbrewer
cols <- brewer.pal(10, 'Spectral')


mean(vals[which(vals<=100)]) # average return period in the 100 year flood zone (21 years)

mean(vals[which(vals>100)]) # average return period in the 500 year flood zone (253 years)


pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse10pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_10pct_interp, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(450, 650, 20)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse2pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_2pct_interp, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(450, 650, 20)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse1pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_1pct_interp, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(450, 650, 20)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_wse02pct.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(wse_02pct_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(450, 650, 20)))
dev.off()
pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_elev.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(elev_agg, xlab = "Easting", ylab = "Northing",col.regions = cols,
             par.settings = list(fontsize = list(text = 20)), at = seq(400, 650, 25)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_probmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(rt_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c("red", "white", "blue")),
             par.settings = list(fontsize = list(text = 20))))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_probmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(prob_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c("red", "white", "blue")),
             par.settings = list(fontsize = list(text = 20)), at = c(10,9,8,7,6,4,5,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_FEMA_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(FEMA_floodbounds, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c("red", "white", "blue")),
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 50)))
dev.off()

pdf("Figures/Muncy_figures/PCHES_Muncy_flood_FEMA_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(FEMA_floodbounds_prob, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c("red", "white", "blue")),
             par.settings = list(fontsize = list(text = 20)), at = c(10,9,8,7,6,4,5,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2)))
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

########### More plots #############

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
pdf("Figures/Selinsgrove_figures/PCHES_zone_vs_interp_boxplot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(Interpolated_RTs~FEMA_zones, data = flood_dat, col = 'lightblue',
        xlab = "FEMA Flood Zones (Return Period)", ylab = "Interpolated Return Period", pch = NA)
dev.off()

flood_prob <- matrix(data = NA, nrow = length(bounds[-allnas]), ncol = 2)
flood_prob[,1] <- 1/bounds[-allnas]*100
flood_prob[,2] <- 1/vals[-allnas]*100
colnames(flood_prob)<-as.character(c("FEMA_prob", "Interpolated_prob"))

# plot FEMA flood zone return periods vs FEMA calculated and interpolated return periods
pdf("Figures/Muncy_figures/PCHES_zone_vs_interp_boxplot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(Interpolated_prob~FEMA_prob, data = flood_prob, col = 'lightblue',
        xlab = "FEMA flood zones AEP", ylab = "Interpolated AEP", pch = NA)
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

############ Goodness of fit analysis and plots ########################
# 1: of all overlapping cells with flood probabilities for each scheme
# what's the average flood probability?

# FEMA flood map
print(c(round(mean(flood_dat[,1])), 'Average FEMA flood zone return period')) # FEMA flood zones: 123 yrs
print(c(round(mean(flood_dat[,2])), 'Average Interpolated flood return period')) # Interpolated flood return periods: 58 yrs

pdf("Figures/Muncy_figures/PCHES_Flood_Interpolation_Histogram.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
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

########################### Root Mean Square Error and more plots ###################################

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
