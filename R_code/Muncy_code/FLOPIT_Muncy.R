###################################################
# file: FLOPIT_Muncy.R
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

rm(list = ls())
graphics.off()

# Install packages
install.packages('raster')
install.packages('rgdal')
install.packages('sp')
install.packages('FNN')
install.packages('dismo')
install.packages('deldir')
install.packages('rgeos')
install.packages('RColorBrewer')
install.packages('scales')
install.packages('gstat')

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

########################## Links to analysis data

# Lidar for Muncy, PA
# http://maps.psiee.psu.edu/ImageryNavigator/

# FEMA flood depth data
# https://msc.fema.gov/portal/advanceSearch

# FEMA log-linear interpolation scheme document
# https://www.fema.gov/media-library-data/1523562952942-4c54fdae20779bb004857f1915236e6c/Flood_Depth_and_Analysis_Grids_Guidance_Feb_2018.pdf

start <- Sys.time()
########################## set working directory
setwd("~/GIS_Data/FLOPIT-master")

source('R_code/FLOPIT_function.R')
########################## Set data links and values
# vector of flood depth raster file locations for Muncy, PA, preprocessed to the same coordinate system and spatial extent
flood_rasters <- c('Data/Muncy_data/depth_10pct_projected_clipped.tif', 'Data/Muncy_data/depth_02pct_projected_clipped.tif',
                   'Data/Muncy_data/depth_01pct_projected_clipped.tif', 'Data/Muncy_data/depth_0_2pct_projected_clipped.tif')
# file location of the elevation map for Muncy, PA, preprocessed to the same coordinate system and spatial extent
elevation_raster <- c('Data/Muncy_data/Muncy_dem_resampled.tif')
# vector of flood probabilities associated with the raster flood depth maps
flood_probabilities <- c(0.1, 0.02, 0.01, 0.002)

######################### FLOPIT analysis: spline interpolation scheme
flopit_Muncy_spline <- FLOPIT(flood_rasters, flood_probabilities, elevation_raster, depth = TRUE, aggregation_value = 1, 
                  method = 'spline', save = FALSE, map_type = 'return period')

######################### FLOPIT analysis: FEMA log-linear interpolation scheme
flopit_Muncy_loglinear <- FLOPIT(flood_rasters, flood_probabilities, elevation_raster, depth = TRUE, aggregation_value = 1, 
                  method = 'log-linear', save = FALSE, map_type = 'return period')

######################## Plots: flood probability and flood zone maps
flopit_map <- flopit_Muncy_spline[[1]] # FLOPIT spline interpolation flood probability map
flopit_zones <- flopit_Muncy_spline[[2]] # flood zone map
end <- Sys.time()

Muncy_time <- end - start

pdf("Figures/Muncy_figures/FLOPIT_Muncy_rtmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 1)))
dev.off()

pdf("Figures/Muncy_figures/FLOPIT_Muncy_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_zones, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
       par.settings = list(fontsize = list(text = 20)), at = seq(0, 500, 1)))
dev.off()

######################## Plots: difference between flood probability map and flood zones
flopit_map_error <- flopit_map
flopit_map_error <- setValues(flopit_map_error, getValues(flopit_zones)-getValues(flopit_map))

pdf("Figures/Muncy_figures/FLOPIT_prob_zone_difference.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_map_error, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
             par.settings = list(fontsize = list(text = 20)), at = seq(-250, 250, 1)))
dev.off()

######################## More plots: flood probability distribution within flood zones

zones <- getValues(flopit_zones)
probabilities <- getValues(flopit_map)

# determine which cells posses NA values in any of the flood return period maps
allnas <- c(which(is.nan(zones)==TRUE), which(is.nan(probabilities)==TRUE),
            which(is.na(zones)==TRUE), which(is.na(probabilities)==TRUE))
allnas <- unique.numeric_version(allnas)

# place cell values with real values across all three return period maps into matrix
flood_dat <- matrix(data = NA, nrow = length(zones[-allnas]), ncol = 2)
flood_dat[,1] <- zones[-allnas]
flood_dat[,2] <- probabilities[-allnas]
colnames(flood_dat)<-as.character(c("FEMA_zones", "Flopit_probabilities"))

# plot FEMA flood zone return periods vs FEMA calculated and interpolated return periods
pdf("Figures/Muncy_figures/zone_vs_probabilities_boxplot.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
boxplot(Flopit_probabilities~FEMA_zones, data = flood_dat, col = 'lightblue',
        xlab = "FEMA Flood Zones (Return Period)", ylab = "Interpolated Return Period", pch = 16, outline = FALSE)
dev.off()