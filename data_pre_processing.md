# Data pre-processing

## Houston

1. Flood water surface elevation and percent annual flood probability raster for Harris county's Sims Bayou watershed were downloaded from the FEMA Flood Map Service Center.
    * FEMA Flood Map Service Center: [https://msc.fema.gov/portal/advanceSearch](https://msc.fema.gov/portal/advanceSearch)
    * *The 4 percent (25 year) flood water surface elevation raster was not used due to elevation quality concerns.*
1. Digital elevation model (DEM) was downloaded for the area at and surrounding Muncy.
    * NED dataset 1/9 arcsecond Digital Elevation Model of Houston Texas: [https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View](https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View)
1. The DEM was reprojected into the coordinate system (TX state plane) and resolution of the FEMA rasters using ArcMap.
1. Flood water surface elevation rasters, percent annual flood probability raster, and dem were clipped to the study area extent.
1. Flood water surface elevation rasters, percent annual flood probability raster, and dem were loaded and processed by the R code.

## Muncy

1. Flood Depth rasters for Lycoming county were downloaded from the FEMA Flood Map Service Center.
    * FEMA Flood Map Service Center: [https://msc.fema.gov/portal/advanceSearch](https://msc.fema.gov/portal/advanceSearch)
1. Digital elevation model (DEM) used in the FEMA flood mapping was downloaded for the area at and surrounding Muncy.
    * PAMAP Program 3.2 ft Digital Elevation Model of Pennsylvania: [http://www.pasda.psu.edu](http://www.pasda.psu.edu)
1. Flood depth rasters were reprojected into the coordinate system (PA state plane) and resolution of the DEM using ArcMap.
1. Flood depth rasters and dem were clipped to the study area extent.
1. Flood depth rasters and dem were loaded and processed by the R code.
