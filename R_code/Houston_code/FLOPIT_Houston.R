###################################################
# file: FLOPIT_Houston.R
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
###################################################

rm(list = ls())
graphics.off()

# Install packages
# install.packages('raster')
# install.packages('rgdal')
# install.packages('sp')
# install.packages('FNN')
# install.packages('dismo')
# install.packages('deldir')
# install.packages('rgeos')
# install.packages('RColorBrewer')
# install.packages('scales')
# install.packages('gstat')

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

# Lidar for Sim bayou, Houston, TX
# https://viewer.nationalmap.gov/basic/

# FEMA flood depth data
# https://msc.fema.gov/portal/advanceSearch

# FEMA log-linear interpolation scheme document
# https://www.fema.gov/media-library-data/1523562952942-4c54fdae20779bb004857f1915236e6c/Flood_Depth_and_Analysis_Grids_Guidance_Feb_2018.pdf

########################## set working directory
source('R_code/FLOPIT_function.R')

########################## Set data links and values
# vector of flood depth raster file locations for Muncy, PA, preprocessed to the same coordinate system and spatial extent
flood_rasters <- c('Input_data/Houston_data/wse_10pct_clipped_resampled.tif', 
                   'Input_data/Houston_data/wse_2pct_clipped_resampled.tif',
                   'Input_data/Houston_data/wse_1pct_clipped_resampled.tif', 
                   'Input_data/Houston_data/wse_02pct_clipped_resampled.tif')
# file location of the elevation map for Muncy, PA, preprocessed to the same coordinate system and spatial extent
elevation_raster <- c('Input_data/Houston_data/ned19_n29x75_w095x50_tx_houstoncity_2008_clipped_feet.tif')
# vector of flood probabilities associated with the raster flood depth maps
flood_probabilities <- c(0.1, 0.02, 0.01, 0.002)

######################### Load FEMA Sims Bayou, Houston, TX flood probability map
str_name <- 'Input_data/Houston_data/pctannchance_clipped_resampled.tif'
prob_map_orig <- raster(str_name)
prob_map_orig <- setValues(prob_map_orig, 1/getValues(prob_map_orig))

######################### FLOPIT analysis: spline interpolation scheme
flopit_Houston_spline <- FLOPIT(flood_rasters, 
                                flood_probabilities, 
                                elevation_raster, 
                                depth = FALSE, 
                                aggregation_value = 1, 
                                method = 'spline', 
                                save_outputs = TRUE, 
                                save_path_map = 'Output/Houston/Data/saved_map_spline',
                                save_path_data = 'Output/Houston/Data/saved_data_spline', 
                                save_path_zones = 'Output/Houston/Data/saved_zones_spline',
                                map_type = 'return period')

######################### FLOPIT analysis: FEMA log-linear interpolation scheme
flopit_Houston_loglinear <- FLOPIT(flood_rasters, 
                                   flood_probabilities, 
                                   elevation_raster, 
                                   depth = FALSE, 
                                   aggregation_value = 1, 
                                   method = 'log-linear', 
                                   save_outputs = TRUE, 
                                   save_path_map = 'Output/Houston/Data/saved_map_logLin',
                                   save_path_data = 'Output/Houston/Data/saved_data_logLin', 
                                   save_path_zones = 'Output/Houston/Data/saved_zones_logLin',
                                   map_type = 'return period')

######################## Plots: flood probability and flood zone maps
flopit_map <- flopit_Houston_spline[[1]] # FLOPIT spline interpolation flood probability map
flopit_zones <- flopit_Houston_spline[[2]] # flood zone map


pdf("Output/Houston/Figures/FLOPIT_Houston_rtmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_map,
             main='FLOPIT Return Period Interpolated Map', 
             xlab = "Easting", 
             ylab = "Northing",
             col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
             par.settings = list(fontsize = list(text = 20)), 
             at = seq(0, 500, 1)))
dev.off()

pdf("Output/Houston/Figures/FLOPIT_Houston_zones.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_zones,
             main='FLOPIT Zones', 
             xlab = "Easting", 
             ylab = "Northing",
             col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
             par.settings = list(fontsize = list(text = 20)), 
             at = seq(0, 500, 1)))
dev.off()

pdf("Output/Houston/Figures/FEMA_probability_map.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(prob_map_orig,
             main='FEMA Probablity Map', 
             xlab = "Easting", 
             ylab = "Northing",
             col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
             par.settings = list(fontsize = list(text = 20)), 
             at = seq(0, 500, 1)))
dev.off()
######################## Plots: difference between flood probability map and flood zones
flopit_map_error <- flopit_map
flopit_map_error <- setValues(flopit_map_error, getValues(flopit_zones)-getValues(flopit_map))

pdf("Output/Houston/Figures/FEMA_FLOPIT_difference.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(flopit_map_error,
             main='Difference Between Interpolated \nand Zonal Return Periods', 
             xlab = "Easting", 
             ylab = "Northing",
             col.regions = colorRampPalette(c('darkred', 'pink', 'skyblue', "blue")),
             par.settings = list(fontsize = list(text = 20)), 
             at = seq(-500, 500, 1)))
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
pdf("Output/Houston/Figures/zone_vs_probabilities_boxplot.pdf",width = 187*0.0393701, height = 187*0.0393701*(100/141))
boxplot(Flopit_probabilities~FEMA_zones, 
        data = flood_dat, 
        col = 'gray',
        horizontal = TRUE, 
        xlim=c(0,3),
        ylab = "FEMA Return Period [years]", 
        xlab = "FLOPIT Interpolated Return Period [years]", pch = 16, outline = FALSE)
dev.off()
