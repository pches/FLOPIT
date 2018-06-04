# PCHES-Flood-Interpolation

This repository is for the PCHES FLOPIT (FLOod Probability Interpolation Tool) project. This tool is distributed under the GNU Public License. No claim of a warranty is expressed or implied, nor does the distribution constitute a warranty.

## Introduction

This document contains code for interpolating inundation return periods for parcels of land between the 10 year and 500 year FEMA flood water surface elevations. The code's goal is to improve flood risk communication and understanding by interpolating flood probabilities from existing FEMA flood maps and data. This analysis focuses on two locations: a model testing case in the Sims Bayou, Houston, TX. and an application case at the town of Muncy, PA.

## Analysis Overview

The code assimilates flood surface elevation or flood depth rasters downloaded from the FEMA flood map services database and digital elevation models, preprocessed in ArcMap to the spatial extent, resolution, and coordinate system necessary for the analysis. It then extrapolates flood surfaces beyond spatial flooding extent, and interpolates inundation return periods for land surface elevations between two flood surface elevations. This analysis produces a raster of interpolated flood return periods over the spatial extent of the study area. All return periods more frequent than the lowest return period are rounded up. All return periods less frequent than the highest return period are determined to be beyond extrapolation range, resulting in NA values.

## Viability at Potential Locations

FLOPIT cannot provide interpolations for areas without FEMA flood water surface elevation or flood depth rasters. These data are not  published on the FEMA flood map service center for every location. As such, FLOPIT is currently restricted by data availability. The following is a non-exhaustive list of locations and FLOPIT's viability at each.

The following cities have the necessary data for FLOPIT.
Houston, TX
Muncy, PA
FLOPIT works in both locations.

FEMA has not published the necessary flood water surface elevation or flood water depth rasters for a FLOPIT analysis for the following cities.
Sellinsgrove, PA
Lewisburg, PA
Tampa, FL
St. Petersburg, FL
Washington, D.C.
San Francisco, CA
Boston, MA
Baltimore, MD
New York, NY

FLOPIT currently (June 4th, 2018) requires either flood water surface elevation or flood water depth rasters of varying return period floods, typically the 500 year, 100 year, 50 year, and 10 year floods. FEMA publishes these data for some locations but not all.

The National Flood Hazard Layer (NFHL), which is published for all the aforementioned cities, includes the spatial extent of the 100 year and 500 year flood zones. If one were to extract the terrain elevation at the edges of each flood zone, one could potentially extrapolate water surface elevations for the 100 year and 500 year floods. These flood elevations, combined with a mean water surface elevation, could be used to interpolation flood probabilities in FLOPIT. This approach is fraut with extrapolation difficulties and is a work in progress.

## Run Instructions

To reproduce this analysis:
1) Download the FLOPIT folder from GitHub
2) In the R_code folder run the Flood_Map_Interpolation_Clean_Muncy.R and Flood_Map_Interpolation_Clean_Houston.R scripts
3) See the outputs in the Figures folder

## Key Plots

![Alt text](/Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_probmap.png)

*In progress... captions and more plots coming soon!*

## Data Overview

**FEMA flood data:** [https://msc.fema.gov/portal/advanceSearch](https://msc.fema.gov/portal/advanceSearch)  
* Muncy data accessed February 8th, 2018
* Houston data accessed January 31st, 2018

**Pennsylvania elevation data:** [http://www.pasda.psu.edu](http://www.pasda.psu.edu)  
* Accessed February 8th, 2018

**Texas elevation data:** [https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View](https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View)  
* Accessed January 31st, 2018

**Texas elevation data used in the FEMA maps:** Discovery Report: Sims Bayou Watershed, HUC 12040104  
* Published June 29th, 2015
* Accessed October 17st, 2017

**2001 Lidar:** [https://coast.noaa.gov/htdata/lidar1_z/geoid12a/data/102/2001_TX_Harris_metadata.html](https://coast.noaa.gov/htdata/lidar1_z/geoid12a/data/102/2001_TX_Harris_metadata.html)

**Supplemented by 2008 Lidar:** [https://tnris.org/data-catalog/entry/houston-galveston-area-council-h-gac-2008-lidar/](https://tnris.org/data-catalog/entry/houston-galveston-area-council-h-gac-2008-lidar/)
* Accessed February 16th, 2018 (data not used in analysis presented here, plans to update analysis with this data are ongoing)

---

*At time of access, all data was free to access (No licenses or paywalls).*
