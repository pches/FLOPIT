###################################################
# file: Selinsgrove_dmg_insurance_mapping.R
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

#rm(list = ls())
#graphics.off()

########################## load packages
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

# Lidar for Selinsgrove, PA
# http://maps.psiee.psu.edu/ImageryNavigator/

# FEMA flood depth data
# https://msc.fema.gov/portal/advanceSearch

# FEMA log-linear interpolation scheme document
# https://www.fema.gov/media-library-data/1523562952942-4c54fdae20779bb004857f1915236e6c/Flood_Depth_and_Analysis_Grids_Guidance_Feb_2018.pdf

########################## set working directory
setwd("~/GIS_Data/FLOPIT-master")

source('R_code/FLOPIT_function_EAD.R') # expected annual damages calculator

source('R_code/FLOPIT_function_NFIP_premiums.R') # estimated annual NFIP premiums calculator

########################## Define input file names and variables

# vector of flood depth raster file locations for Selinsgrove, PA, preprocessed to the same coordinate system and spatial extent
flood_rasters <- c('Data/Sellinsgrove_data/dg_10pct_clip.tif', 'Data/Sellinsgrove_data/dg_4pct_clip.tif',
                   'Data/Sellinsgrove_data/dg_2pct_clip.tif', 'Data/Sellinsgrove_data/dg_1pct_clip.tif',
                   'Data/Sellinsgrove_data/dg_02pct_clip.tif')
# file location of the elevation map for Selinsgrove, PA, preprocessed to the same coordinate system and spatial extent
elevation_raster <- c('Data/Sellinsgrove_data/Sellinsgrove_DEM_clip1.tif')
# vector of flood probabilities associated with the raster flood depth maps
flood_probabilities <- c(0.1, 0.04, 0.02, 0.01, 0.002)


damage_val <- c(0.20,0.44,0.58,0.68,0.78,0.85,0.92,0.96,1,1,1,1,1,1,1,1,1,1)
depth_val <- c(0,0.5,1,1.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15)*3.28084


pointstr <- 'Data/Sellinsgrove_data/Selinsgrove_properties_project1.shp'
selinsgrove_buildings <- shapefile(pointstr)

home_elev <- mask(x = raster(elevation_raster), mask = selinsgrove_buildings, maskvalue=NA)

home_elev <- extract(x = raster(elevation_raster), y = selinsgrove_buildings)

elev <- raster(elevation_raster)

coordinates(selinsgrove_buildings)

building_grid <- selinsgrove_buildings

gridded(building_grid)=TRUE

building_grid <- rasterToPoints(selinsgrove_buildings, spatial = TRUE)

########################## Expected annual damages as percent of home value
ann_dmg <- FLOPIT_EAD(flood_rasters, flood_probabilities, elevation_raster, depth = TRUE, aggregation_value = 10, 
                   method = 'log-linear', map_type = 'return period', save = FALSE, 
                   damage_vals = damage_val, depth_vals = depth_val)

########################## Estimated NFIP premiums as percent of home value
NFIP_insurance <- FLOPIT_NFIP_premiums(flood_rasters, flood_probabilities, elevation_raster, depth = TRUE, aggregation_value = 10, 
                   method = 'log-linear', map_type = 'return period', save = FALSE, Coverage_Building = 100000,
                   Coverage_Contents = 60000)

########################## Figures
expected_damages_map <- ann_dmg[[3]]
insurance_rate_map <- NFIP_insurance[[3]]

building_val_assess <- as.integer(selinsgrove_buildings$assdttlval)


building_dmg <- extract(x = expected_damages_map, y = selinsgrove_buildings)

building_NFIP <- extract(x = insurance_rate_map, y = selinsgrove_buildings)


pdf("Figures/Selinsgrove_figures/FLOPIT_Selinsgrove_EAD.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(expected_damages_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('blue', 'skyblue', 'pink', 'darkred')),
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 16, 0.1)))
dev.off()

pdf("Figures/Selinsgrove_figures/FLOPIT_Selinsgrove_NFIP.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(insurance_rate_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('blue', 'skyblue', 'pink', 'darkred')),
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 16, 0.1)))
dev.off()

pdf("Figures/Selinsgrove_figures/FLOPIT_Selinsgrove_DEM.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(elev, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(terrain.colors(100)),
             par.settings = list(fontsize = list(text = 20)), at = seq(400, 550, 0.1)))
dev.off()

diff_map <- insurance_rate_map
diff_map <- setValues(diff_map, (getValues(insurance_rate_map)-getValues(expected_damages_map)))

pdf("Figures/Selinsgrove_figures/FLOPIT_Selinsgrove_diffmap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
print(spplot(diff_map, xlab = "Easting", ylab = "Northing",col.regions = colorRampPalette(c('blue', 'skyblue', 'pink', 'darkred')),
             par.settings = list(fontsize = list(text = 20)), at = seq(0, 16, 0.1)))
dev.off()

###################

pdf("Figures/Selinsgrove_figures/FLOPIT_Selinsgrove_building_map.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
plot(selinsgrove_buildings, pch = 15, cex = 0.7, 
     col = colorRampPalette(c('blue', 'skyblue', 'pink', 'darkred'))(12)[(building_NFIP-building_dmg)+1]
     #,add=TRUE
     )
legend('topleft', 
        c("0","1","2",'3','4','5','6','7','8','9','10','11'), fill = colorRampPalette(c('blue', 'skyblue', 'pink', 'darkred'))(12), xpd = NA)
dev.off()

#building_dmg[which(building_dmg==0)]<-0.02

selinsgrove_buildings$risk = building_dmg
selinsgrove_buildings$NFIP = building_NFIP
selinsgrove_buildings$r_N = building_NFIP-building_dmg

building_dmg[which(building_dmg==0)]<- NA
selinsgrove_buildings$percent = building_NFIP/building_dmg*100

  
building_NFIP[which(is.na(building_dmg)==TRUE)]<-NA

total_actuarial <- sum(building_dmg/100*building_val_assess, na.rm = TRUE)

total_NFIP <- sum(building_NFIP/100*building_val_assess, na.rm = TRUE)

1/(total_actuarial/total_NFIP)*100

writeOGR(selinsgrove_buildings, "Data/Sellinsgrove_data", "selinsgrove_buildings_withdmg7", driver="ESRI Shapefile")
