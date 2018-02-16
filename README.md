# PCHES-Flood-Interpolation
This repository is for the PCHES FLOPIT (FLOod Probability Interpolation Tool) project.
This tool is distributed under the GNU Public License.
No claim of a warranty is expressed or implied, nor does the distribution constitute a warranty.

Introduction
This document contains code for interpolating flood return periods between the 10 year and 500 year FEMA floods. 
The code's goal is to improve flood risk communication and understanding by interpolating flood probabilities from existing FEMA flood maps and data. This analysis focuses on two locations: model testing case in the Sims Bayou, Houston, TX. and an application case at the town of Muncy, PA.

Analysis Overview
The code assimilates flood surface elevation or flood depth rasters downloaded from the FEMA flood map services database and digital elevation models, preprocessed in ArcMap to the spatial extent, resolution, and coordinate system necessary for the analysis. It then extrapolates flood surfaces beyond spatial flooding extent, and interpolates flood return periods for land surface elevations between two flood surface elevations. This analysis produces a raster of interpolated flood return periods over the spatial extent of the study area. All return periods more frequent than the lowest return period are rounded up. All return periods less frequent than the highest return period are determined to be beyond extrapolation range, resulting in NA values.

Key Plots

(https://github.com/Joelroopeckart/PCHES-Flood-Interpolation/Figures/Muncy_figures/PCHES_Muncy_flood_interpolation_probmap.png)


Data Overview
FEMA flood data: https://msc.fema.gov/portal/advanceSearch
Muncy data accessed February 8th, 2018
Houston data accessed January 31st, 2018

Pennsylvania elevation data: http://www.pasda.psu.edu
Accessed February 8th, 2018

Texas elevation data: https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View
Accessed January 31st, 2018

Texas elevation data used in the FEMA maps: 
Discovery Report: Sims Bayou Watershed, HUC 12040104
Published June 29th, 2015
Accessed October 17st, 2017
2001 Lidar https://coast.noaa.gov/htdata/lidar1_z/geoid12a/data/102/2001_TX_Harris_metadata.html
Suplimented by 2008 Lidar: https://tnris.org/data-catalog/entry/houston-galveston-area-council-h-gac-2008-lidar/

Accessed February 16th, 2018 (data not used in analysis presented here, plans to update analysis with this data are ongoing)

At time of access, all data was free to access (No licenses or paywalls).
