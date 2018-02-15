# PCHES-Flood-Interpolation
This repository is for the PCHES flood map interpolation project

Introduction
This document contains code for interpolating flood return periods between the 10 year and 500 year FEMA floods. 
The code's goal is to improve flood risk communication and understanding by interpolating flood probabilities from existing FEMA flood maps and data.
This analysis focuses on two locations: model testing case in the Sims Bayou, Houston, TX. and an application case at the town of Muncy, PA.

Analysis Overview
The code assimilates flood surface elevation or flood depth rasters downloaded from the FEMA flood map services database and digital elevation models, preprocessed in ArcMap to the spatial extent, resolution, and coordinate system necessary for the analysis. It then extrapolates flood surfaces beyond spatial flooding extent, and interpolates flood return periods for land surface elevations between two flood surface elevations. This analysis produces a raster of interpolated flood return periods over the spatial extent of the study area. All return periods more frequent than the lowest return period are rounded up. All return periods less frequent than the highest return period are determined to be beyond extrapolation range, resulting in NA values.

FEMA flood data: https://msc.fema.gov/portal/advanceSearch

Pennsylvania elevation data: http://www.pasda.psu.edu

Texas elevation data: https://viewer.nationalmap.gov/basic/?basemap=b1&category=ned,nedsrc&title=3DEP%20View
