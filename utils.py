'''
This code is the code that was used for the analysis of this research paper:

The FLOod Probability Interpolation Tool (FLOPIT): A Simple
Tool to Improve Spatial Flood Probability Quantification
and Communication

Authors: Mahkameh Zarekarizi (1,2,*), K. Joel Roop-Eckart (3) , Sanjib Sharma (1), and Klaus Keller (1,3)

1: Earth and Environmental Systems Institute, Pennsylvania State University, University Park, PA 16802, USA
2: Jupiter Intelligence, San Mateo, CA 94401, USA
3: Department of Geosciences, Pennsylvania State University, University Park, PA 16802, USA
*Author to whom correspondence should be addressed (mahkameh.zare@gmail.com).

Author and copyright: Mahkameh Zarekarizi, Jupiter Intelligence, and Pennsylvania State University  
This notebook is written by Mahkameh Zarekarizi, while at Jupiter Intelligence 

Distributed under the GNU general public license

Disclaimer: The datasets, software tools, results, and any other resources associated with FLOPIT
and this manuscript are intended for academic research and education (not for real-world decisionmaking) and provided as-is without warranty of any kind, express or implied. In no event shall the
authors or copyright holders be liable for any claim, damages, or other liability in connection with
the use of these resources

Acknowledgement: Special thanks to Dr. Luke Madaus for contributing to this code by writing the function for filling the NaNs with their nearest neighbor that is not NaN
'''
    
# Import
import xarray
import rasterio
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import KDTree
from rasterio import crs
import os 
from scipy import interpolate
from sklearn.metrics import auc

def log_linear_interpolate(x, y, xout, kind='linear',plot=False,verbose=1):
    
    if (np.sum(np.isnan(y)) == len(y)) or (np.sum(y<=0) == len(y)):
        yout = np.nan
        return yout
    
    if (np.sum(x<=0) > 1) or (np.sum(y<=0) > 1) or (xout<=0):
        # do a simple linear regression 
        # In coastal watersheds, elevation could become negative
        if plot:
            print('NOTE: Switching to a simple regression because elevation is below 0') if verbose == 1 else None
        interp_function = scipy.interpolate.interp1d(x, y, kind=kind,assume_sorted=True,bounds_error=False,fill_value="extrapolate")
        yout = interp_function(xout)
        
        if plot:
            plt.figure()
            plt.plot(x,y, 'o-',color='blue')
            plt.scatter([xout],[yout],color='red')
            print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
        
        
    else:
        logx = np.log10(x)
        #logy = np.log10(y)
        interp_function = scipy.interpolate.interp1d(logx, y, kind=kind,assume_sorted=True,bounds_error=False,fill_value="extrapolate")
        yout = interp_function(np.log10(xout))
        
        if plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.plot(logx,y, 'o-',color='blue')
            plt.scatter([np.log10(xout)],[interp_function(np.log10(xout))],color='red')
            
            plt.subplot(1,2,2)
            plt.plot(x,y, 'o-',color='blue')
            plt.scatter([xout],[yout],color='red')
            
            print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
    
    return yout

def linear_interpolate(x, y, xout, kind='linear',plot=False,verbose=1):

    if (np.sum(np.isnan(y)) == len(y)) or (np.sum(y<=0) == len(y)):
        yout = np.nan
        return yout

    interp_function = scipy.interpolate.interp1d(x, y, kind=kind,assume_sorted=True,bounds_error=False,fill_value="extrapolate")
    yout = interp_function(xout)
        
    if plot:
        plt.figure()
        plt.plot(x,y, 'o-',color='blue')
        plt.scatter([xout],[yout],color='red')
        print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
    
    return yout


def spline_interpolate(x,y,xout,plot=False,verbose=1):
    
    if (np.sum(np.isnan(y)) == len(y)) or (np.sum(y<=0) == len(y)):
        yout = np.nan
        return yout
    
    increasing = all(i < j for i, j in zip(x, x[1:]))
    if increasing:
        cs = CubicSpline(x, y)
        yout = cs(xout)
        
        if plot:
            plt.figure()
            plt.plot(x,y, 'o-',color='blue')
            plt.scatter([xout],[yout],color='red')
            xplot = np.linspace(np.min(x),np.max(x),num=100)
            plt.plot(xplot,cs(xplot),'.',color='yellow')
            print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
        
    else:
        yout = None 
    
    
    return yout


def log_spline_interpolate(x, y, xout, kind='linear',plot=False,verbose=1):
        
    # If all are NaN or all are negative, return NaN
    if (np.sum(np.isnan(y)) == len(y)) or (np.sum(y<=0) == len(y)):
        yout = np.nan
        if plot:
            print('For this sample point all data are NaN or negative. Returning NaN')
        return yout
    
    # If there are negatives in x, y, or xout, switch to simple regression
    if (np.sum(x<=0) > 1) or (np.sum(y<=0) > 1) or (xout<=0):
        
        if plot:
            print('NOTE: We are switching to a simple regression because elevation is below 0') if verbose == 1 else None if verbose == 1 else None
        
        interp_function = scipy.interpolate.interp1d(x, y, kind=kind,assume_sorted=True,bounds_error=False,fill_value="extrapolate")
        yout = interp_function(xout)
        
        if plot:
            plt.figure()
            plt.plot(x,y, 'o-',color='blue')
            plt.scatter([xout],[yout],color='red')
            print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
        
        
    else:
        logx = np.log10(x)
        logy = np.log10(y)
        increasing = all(i < j for i, j in zip(logx, logx[1:])) 
        if increasing and (np.all(np.isfinite(logy))):
            interp_function = CubicSpline(logx, logy)
            yout = np.power(10.0, interp_function(np.log10(xout)))
        
            if plot:
                plt.figure()
                plt.subplot(1,2,1)
                plt.plot(logx,logy, 'o-',color='blue')
                plt.scatter([np.log10(xout)],[interp_function(np.log10(xout))],color='red')
            
                plt.subplot(1,2,2)
                plt.plot(x,y, 'o-',color='blue')
                plt.scatter([xout],[yout],color='red')
                
                print(f'Data received for the sample point: x={x}, y={y}, xout={xout}, yout={yout}') if verbose == 1 else None
        else: 
            if plot:
                print('log(x) is not increasing or there is an +/-Inf in log(y). Returning NaN')

            yout = np.nan
    
    return yout

def generate_tif(src,
    output: np.ndarray,
    outfile_name: str,
    ) -> None:
    """
    Creates a tiff file with the given Affine.

    :param use_affine: The affine to use.
    :param output: The data to write to the tiff. Must be 2D or 3D.
    :param outfile_name: The resulting filename.
    :return:
    """
    
    # get info based on shape of array
    width = output.shape[1]
    height = output.shape[0]
    count = 1
    output = np.array([output]) 

    with rasterio.open(src, "r") as src:
        input_transform = src.transform

    profile = {
        "driver": "GTiff",
        "dtype": rasterio.dtypes.float32,
        "nodata": -9999,
        "width": width,
        "height": height,
        "count": count,
        "transform": input_transform,
    }
    with rasterio.open(outfile_name, "w", **profile) as dst:
        for band in range(output.shape[0]):
            dst.write(output[band].astype(rasterio.dtypes.float32), band + 1)
    return

def plot_dem(elev,input_dem_plot_path,units='m',sample_point=None):
    plt.figure(figsize=(12,8))
    plt.pcolor(elev.x,elev.y,elev)
    plt.colorbar()
    if sample_point is not None:
        plt.scatter(sample_point[0],sample_point[1],color='red')
    plt.title(f'DEM({units})')
    plt.savefig(input_dem_plot_path)
    plt.close()

def plot_flood_depths(raster_data_ds,input_depth_plot_path,units='m',type='depth',sample_point=None):
    
    plt.figure(figsize=(12,8))

    for i in range(len(raster_data_ds.returnPeriod)):
        plt.subplot(2,2,i+1)
        plt.pcolor(raster_data_ds.x,raster_data_ds.y,raster_data_ds.isel(returnPeriod=i))
        if sample_point is not None:
            plt.scatter(sample_point[0],sample_point[1],color='red')

        plt.colorbar()
        plt.title(f'water surface elevation with {raster_return_periods[i]}% chance')
    plt.savefig(input_depth_plot_path)
    plt.close()

def generate_depth_elev_ds(flood_rasters_path_list:list, #list of paths to the raster data  
           input_return_periods:list, #list of floats that contains exceedence probabilities of the provided rasters. This must be the same length as the first argument  
           elevation_raster: str, #path to the elevation DEM data 
           depth_or_wse:str, #If the provided data are depth, set this as 'depth'. If they are water surface elevation, this this to 'wse'. 
           depth_units:str, #The units for the input flood rasters. Options are 'ft', 'feet', 'm', 'meter', and 'meters'
           dem_units: str, #The units for the DEM data. Options are 'ft', 'feet', 'm', 'meter', and 'meters'
           verbose:int = 1,
        ):
    # Count the number of input rasters 
    Nrasters = len(flood_rasters_path_list)
    print(f"Received {Nrasters} rasters") if verbose == 1 else None

    if Nrasters != len(input_return_periods):
        print(f"Error: number of return periods must be equal to the number of input rasters") if verbose == 1 else None
    
    raster_exc_probs = 1/np.array(input_return_periods)

    # import Lidar DEM of the study area
    dem_raster = xarray.open_rasterio(elevation_raster)
    elev = dem_raster.isel(band=0)
    del elev.coords['band']

    # Read flood raster data in a list  
    flood_raster_list = []
    for i in range(Nrasters):
        print(f'reading raster {flood_rasters_path_list[i]}') if verbose == 1 else None
        path = flood_rasters_path_list[i]
        tmp_raster = xarray.open_rasterio(path)
        del tmp_raster.coords['x']
        del tmp_raster.coords['y']
        tmp_raster.coords['x'] = elev.x
        tmp_raster.coords['y'] = elev.y
                
        if depth_units[i] in ['ft','feet']:
            print(f'Received data for raster {flood_rasters_path_list[i]} are in "feet". Converting them to meters') if verbose == 1 else None
            tmp_raster = tmp_raster * 0.3048
        elif depth_units[i] in ['meters','m','meter']:
            print(f'Received data for raster {flood_rasters_path_list[i]} are in "meters".') if verbose == 1 else None
        else:
            print(f'Error: Input units for {flood_rasters_path_list[i]} is not known') if verbose == 1 else None
        
        
        flood_raster_list.append(tmp_raster)
        
    # Convert the list to a dataset by concatenating against a new dimension that shows exccedance probability
    raster_data_ds = xarray.concat(flood_raster_list,dim='returnPeriod')
    
    # Getting rid of the "band" dimension 
    raster_data_ds = raster_data_ds.isel(band=0)
    del raster_data_ds.coords['band']
    
    # Make the exceedance probability a new coordinate and assign values
    raster_data_ds.coords['returnPeriod']=input_return_periods
        
    # TODO: make sure the x and y coordinates align well
    del elev.coords['y']
    elev.coords['y']=raster_data_ds.coords['y']

    del elev.coords['x']
    elev.coords['x']=raster_data_ds.coords['x']
            
    if dem_units in ['ft','feet']:
        print('Received DEM data are in "feet". Converting them to meters') if verbose == 1 else None
        elev = elev * 0.3048        
    elif dem_units in ['meters','m','meter']:
        print('Received DEM data are in "meters".') if verbose == 1 else None
    else:
        print('Error: Input units for DEM is not known') if verbose == 1 else None

    # Convert flood depths to elevations if needed
    
    if depth_or_wse == 'depth':
        print('Received data are water depths... Adding elevation of depth values') if verbose == 1 else None
        wse = raster_data_ds.copy(deep=True)
        depth = raster_data_ds
        depth = depth.where(depth>=0,np.nan)
        
        # convert depths to elevations
        for i in range(Nrasters):
            wse[i,:,:] = raster_data_ds[i,:,:] + elev.values
            
    elif depth_or_wse == 'wse':
        print('Received data are water surface elevations...') if verbose == 1 else None
        wse = raster_data_ds
        depth = raster_data_ds.copy(deep=True)
        
        for i in range(Nrasters):
            depth[i,:,:] = raster_data_ds[i,:,:] - elev.values
    else:
        print('Error: argument "depth_or_wse" recived a string that is not in the options list. Options are "depth" or "wse"') if verbose == 1 else None

    # If depth is negative, make them NaN 
    depth = depth.where(depth>=0,np.nan)

    # If depth is negative, make them NaN
    wse = wse.where(wse>=-100000,np.nan)
    
    flood_ds = depth.to_dataset(name='depth')
    flood_ds = flood_ds.assign(wse=wse)
    flood_ds = flood_ds.assign(elev=elev)
    
    return flood_ds

def point_EAD(DD_depth,
              DD_damage,
              flood_rasters_path_list:list,
              elevation_raster_path:str,
              raster_return_periods:list,
              depth_or_wse:str,
              point_lat:float,
              point_lon:float,
              output_plot_dir:str,
              Struc_Value:float,
              depth_units:str,
              dem_units:str,
              verbose:int = 1,
             ):
    
    plt.figure()
    plt.plot(DD_depth,DD_damage,'o-')
    plt.xlabel('Depth(m)')
    plt.ylabel('Damage')
    plt.grid()
    plt.savefig(output_plot_dir + 'depth_damage_function.png')
    plt.close()

    ds = generate_depth_elev_ds(flood_rasters_path_list=flood_rasters_path_list,
               input_return_periods = raster_return_periods,
               elevation_raster = elevation_raster_path,
               depth_or_wse = depth_or_wse,
               depth_units = depth_units,
               dem_units = dem_units)
    
    #print(ds.wse.sel(x=3144641.69354839,y=13810613.9516129,method='nearest').values)
    #print(ds.depth.sel(x=3144641.69354839,y=13810613.9516129,method='nearest').values)
    damage_vals = DD_damage * Struc_Value

    input_flood_depth = ds.depth
    point_depths = input_flood_depth.sel(x=point_lon,y=point_lat,method='nearest').values
    print(f'Flood depth at the requested point are:{point_depths} for these return periods:{raster_return_periods}') if verbose == 1 else None
    
    plt.figure()
    plt.plot(raster_return_periods,point_depths,'o-')
    plt.xlabel('Retrn Period(yr)')
    plt.ylabel('WSE(m)')
    plt.grid()
    plt.savefig(output_plot_dir + 'depth_return_period_curve.png')

    plt.figure()
    plt.plot(1/raster_return_periods,point_depths,'o-')
    plt.xlabel('Exc. Prob.')
    plt.ylabel('Depth(m)')
    plt.grid()
    plt.savefig(output_plot_dir + 'depth_exc_prob_curve.png')
    plt.close()


    InterpFunction = interpolate.interp1d(DD_depth,damage_vals)
    point_damages = InterpFunction(point_depths)
    
    plt.figure()
    plt.plot(1/raster_return_periods,point_damages,'o-')
    plt.xlabel('Exc. Prob.')
    plt.ylabel('Damage(USD)')
    plt.grid()
    plt.title('EPL Curve')
    plt.savefig(output_plot_dir + 'damage_exc_prob_curve.png')
    plt.close()

    ead = auc(1/raster_return_periods,point_damages)
    print(f'Annual Expected Damages (EAD)= ${ead}') if verbose == 1 else None
    return ead


################################### FLOPIT function ##########################################
def FLOPIT(flood_rasters_path_list:list, #list of paths to the raster data  
           input_return_periods:list, #list of floats that contains exceedence probabilities of the provided rasters. This must be the same length as the first argument  
           elevation_raster: str, #path to the elevation DEM data 
           depth_or_wse:str, #If the provided data are depth, set this as 'depth'. If they are water surface elevation, this this to 'wse'. 
           depth_units:str, #The units for the input flood rasters. Options are 'ft', 'feet', 'm', 'meter', and 'meters'
           dem_units: str, #The units for the DEM data. Options are 'ft', 'feet', 'm', 'meter', and 'meters'
           interpolation_method:str, #The method for interpolation. Options are 'log-linear' or 'spline' 
           flood_zone_return_period:list, #If you want to output data/maps for a specific return period, enter the return period here. Otherwise, set as 0
           output_data_dir: str,
           output_plot_dir: str,
           aep_raster_generate:int, #The path where you want to save the AEP data. Will be saved as raster (.tif). If you do not want this, set it as None
           aep_map_generate:int, #The path where you want to save the AEP map plots. Will be saved as .png. If you do not want this, set it as None
           flood_zone_raster_generate:int, #If you want to output data/maps for a specific return period, enter the path where you want to save raster data (as .tif) 
           flood_zone_map_generate:int, #If you want to output data/maps for a specific return period, enter the path where you want to save the maps (as .png) 
           extrapolate: int = 1, #If the interpolation should be allowed to do extrapolation. default is 1 (yes). If no, set it to 0
           input_dem_plot_generate: int=0, #If you want to save the plot showing the DEM map, set the path here. Otherwise, set to None
           input_depth_plot_generate: int=0, #If you want to save the plot showing the flood maps, set the path here. Otherwise, set to None
           sample_point = None,
           backfill_nans=False,
           areaName = 'None',
           verbose:str = 1,
           
          ):
    '''
    This function returns the xarray dataset that contains the interpolated probabilities
    '''
    
    if not os.path.isdir(output_plot_dir):
        print(f"creating a directory for output plots: {output_plot_dir}") if verbose == 1 else None
        os.mkdir(output_plot_dir)
    if not os.path.isdir(output_data_dir):
        print(f"creating a directory for output plots: {output_data_dir}") if verbose == 1 else None
        os.mkdir(output_data_dir)
                
    ds = generate_depth_elev_ds(
            flood_rasters_path_list = flood_rasters_path_list,
            input_return_periods = input_return_periods,
            elevation_raster = elevation_raster,
            depth_or_wse = depth_or_wse,
            depth_units = depth_units,
            dem_units = dem_units,
            verbose = verbose)
    
    n_inputs = len(flood_rasters_path_list)
    
    elev = ds.elev
    input_flood_wse = ds.wse
    input_flood_wse_tmp = ds.wse
    
    raster_exc_probs = 1/np.array(input_return_periods)

    if input_dem_plot_generate == 1:
        print(f'Plotting and saving the input DEM map to {output_plot_dir}/input_dem.png') if verbose == 1 else None
        plot_dem(elev, output_plot_dir + '/input_dem.png',units='m',sample_point=sample_point)
        
    if input_depth_plot_generate == 1:
        print(f'Plotting and saving the input flood maps to {output_plot_dir + "input_depths.png"}') if verbose == 1 else None
        plot_flood_depths(input_flood_wse, output_plot_dir + '/input_depths.png',units='m',sample_point=sample_point)
    
    # Fill the NaNs with their nearest neighbor values 
    print('Backfilling the NaNs with data from nearest neighbor cell...') if verbose == 1 else None
    if backfill_nans:
        input_flood_wse_backfilled = backfill_nan_with_nearest(input_flood_wse)
        input_flood_wse_backfilled = input_flood_wse_backfilled.where(ds.wse.isel(returnPeriod=n_inputs-1)>-1000,np.nan)
    else:
        input_flood_wse_backfilled = input_flood_wse 

    # Because the backfilling function fills NaNs in the oroginal flood data as well.  
    input_flood_wse = input_flood_wse_tmp
    
    bc = xarray.broadcast(input_flood_wse_backfilled.returnPeriod,input_flood_wse_backfilled)
    bc_rp = bc[0]
    bc_aep = 1/bc_rp

    # Start the interpolation
    if interpolation_method == 'spline':
        print('Starting the interpolation using the "Spline" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            spline_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = spline_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
        
    elif interpolation_method == 'log-linear':
        print('Starting the interpolation using the "Log-Linear" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            log_linear_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = log_linear_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
    elif interpolation_method == 'linear':
        print('Starting the interpolation using the "Linear" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            linear_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = linear_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)

    elif interpolation_method == 'log-spline':
        print('Starting the interpolation using the "Log-Spline" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            log_spline_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = log_spline_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)

    print('Finished the interpolation') if verbose == 1 else None
    interpolated_aep = interpolated_aep.load()
    interpolated_aep = interpolated_aep.where(input_flood_wse.isel(returnPeriod=n_inputs-1) >= -1000,np.nan)
    
    if extrapolate == 0:
        interpolated_aep = interpolated_aep.where(interpolated_aep <= np.max(raster_exc_probs),np.max(raster_exc_probs))
        interpolated_aep = interpolated_aep.where(interpolated_aep >= np.min(raster_exc_probs),np.min(raster_exc_probs))
        
        print(f'Converted AEP values greater than {np.max(raster_exc_probs)} to {np.max(raster_exc_probs)} and values less than {np.min(raster_exc_probs)} to NaN') if verbose == 1 else None
        
    if aep_raster_generate == 1:
        generate_tif(elevation_raster,interpolated_aep,output_data_dir + f'/aep_{interpolation_method}_{areaName}.tif')
        
    if aep_map_generate == 1:
        print('Saving the AEP plot to ...') if verbose == 1 else None
        plt.figure(figsize=(20,15))
        plt.pcolor(interpolated_aep.x,interpolated_aep.y,interpolated_aep,vmin=np.min(raster_exc_probs),vmax=np.max(raster_exc_probs))
        plt.colorbar()
        plt.title('AEP')
        if sample_point is not None:
            plt.scatter(sample_point[0],sample_point[1],color='red')
        plt.savefig(output_plot_dir + f'/aep_map_{interpolation_method}_{areaName}.png')
        plt.close()
        
    if flood_zone_return_period is not None:
        n_output_rp = len(flood_zone_return_period)
        print(f'Generating output for {n_output_rp} new return periods:{flood_zone_return_period}') if verbose == 1 else None
        print('Checking the validity of the return periods...') if verbose == 1 else None
        
        if (np.min(flood_zone_return_period) < np.min(input_return_periods)) and (extrapolate == 0):
            print('Error: extrapolate is set to zero but a return period outside of the input return period range is requested') if verbose == 1 else None
        elif (np.max(flood_zone_return_period) > np.max(input_return_periods)) and (extrapolate == 0):
            print('Error: extrapolate is set to zero but a return period outside of the input return period range is requested') if verbose == 1 else None

        flood_zone_return_period = np.array(flood_zone_return_period)
        output_exc_prob = 1/flood_zone_return_period
        
        for j in range(0,len(flood_zone_return_period)):
            output_rp_zone = interpolated_aep.where(interpolated_aep >= output_exc_prob[j], np.nan)
            
            if flood_zone_raster_generate == 1:
                generate_tif(elevation_raster,output_rp_zone,output_data_dir + f'/estimated_aep_data_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.tif')
                
            dataC = np.zeros(elev.shape,float)
            dataC[:] = flood_zone_return_period[j]
            return_period_da = elev.copy(deep=True,data=dataC)

            # Start the interpolation
            if interpolation_method == 'spline':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Spline method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            spline_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = spline_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None
                    
            # TODO: add log-linear code 
            elif interpolation_method == 'log-linear':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Log-Linear method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            log_linear_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = log_linear_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

            elif interpolation_method == 'linear':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Linear method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            linear_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = linear_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

            elif interpolation_method == 'log-spline':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Log-Spline method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            log_spline_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = log_spline_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

    
            interpolated_wse = interpolated_wse.where(interpolated_aep >= output_exc_prob[j], np.nan)
            generate_tif(elevation_raster,interpolated_wse.values,output_data_dir + f'/estimated_wse_data_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.tif')
            print(f'Finished the interpolation of water surface elevation for the requested new return period: {flood_zone_return_period[j]}') if verbose == 1 else None

            
            if flood_zone_map_generate == 1:
                plt.figure(figsize=(20,15))
                plt.pcolor(output_rp_zone.x,output_rp_zone.y,output_rp_zone,vmin=np.min(output_rp_zone),vmax=np.max(output_rp_zone))
                plt.colorbar()
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood AEP')
                plt.savefig(output_plot_dir + f'/estimated_aep_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
                
                output_rp_aep_inout = output_rp_zone.where(np.isnan(output_rp_zone),1)
                plt.figure(figsize=(20,15))
                plt.pcolor(output_rp_aep_inout.x,output_rp_aep_inout.y,output_rp_aep_inout,vmin=0,vmax=1/100)
                plt.colorbar()
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood in-out zone')
                plt.savefig(output_plot_dir + f'/estimated_aep_zonal_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
                
                plt.figure(figsize=(20,15))
                plt.pcolor(interpolated_wse.x,interpolated_wse.y,interpolated_wse)
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.colorbar()
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood WSE')
                plt.savefig(output_plot_dir + f'/estimated_wse_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
    return
    
################################### FLOPIT function ##########################################
def FLOPIT_preProcessed(ds, 
           input_return_periods:list, #list of floats that contains exceedence probabilities of the provided rasters. This must be the same length as the first argument  
           elevation_raster: str, #path to the elevation DEM data 
           interpolation_method:str, #The method for interpolation. Options are 'log-linear' or 'spline' 
           flood_zone_return_period:list, #If you want to output data/maps for a specific return period, enter the return period here. Otherwise, set as 0
           output_data_dir: str,
           output_plot_dir: str,
           aep_raster_generate:int, #The path where you want to save the AEP data. Will be saved as raster (.tif). If you do not want this, set it as None
           aep_map_generate:int, #The path where you want to save the AEP map plots. Will be saved as .png. If you do not want this, set it as None
           flood_zone_raster_generate:int, #If you want to output data/maps for a specific return period, enter the path where you want to save raster data (as .tif) 
           flood_zone_map_generate:int, #If you want to output data/maps for a specific return period, enter the path where you want to save the maps (as .png) 
           extrapolate: int = 1, #If the interpolation should be allowed to do extrapolation. default is 1 (yes). If no, set it to 0
           input_dem_plot_generate: int=0, #If you want to save the plot showing the DEM map, set the path here. Otherwise, set to None
           input_depth_plot_generate: int=0, #If you want to save the plot showing the flood maps, set the path here. Otherwise, set to None
           sample_point = None,
           areaName = 'None',
           verbose:str = 1,
           
          ):
    '''
    This function returns the xarray dataset that contains the interpolated probabilities
    '''
    
    if not os.path.isdir(output_plot_dir):
        print(f"creating a directory for output plots: {output_plot_dir}") if verbose == 1 else None
        os.mkdir(output_plot_dir)
    if not os.path.isdir(output_data_dir):
        print(f"creating a directory for output plots: {output_data_dir}") if verbose == 1 else None
        os.mkdir(output_data_dir)
    
    n_inputs = len(input_return_periods)
    
    elev = ds.elev
    input_flood_wse = ds.wse
    input_flood_wse_tmp = ds.wse
    
    raster_exc_probs = 1/np.array(input_return_periods)

    if input_dem_plot_generate == 1:
        print(f'Plotting and saving the input DEM map to {output_plot_dir}/input_dem.png') if verbose == 1 else None
        plot_dem(elev, output_plot_dir + '/input_dem.png',units='m',sample_point=sample_point)
        
    if input_depth_plot_generate == 1:
        print(f'Plotting and saving the input flood maps to {output_plot_dir + "input_depths.png"}') if verbose == 1 else None
        plot_flood_depths(input_flood_wse, output_plot_dir + '/input_depths.png',units='m',sample_point=sample_point)
    
    input_flood_wse_backfilled = input_flood_wse 

    # Because the backfilling function fills NaNs in the oroginal flood data as well.  
    input_flood_wse = input_flood_wse_tmp
    
    bc = xarray.broadcast(input_flood_wse_backfilled.returnPeriod,input_flood_wse_backfilled)
    bc_rp = bc[0]
    bc_aep = 1/bc_rp

    # Start the interpolation
    if interpolation_method == 'spline':
        print('Starting the interpolation using the "Spline" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            spline_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = spline_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
        
    elif interpolation_method == 'log-linear':
        print('Starting the interpolation using the "Log-Linear" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            log_linear_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = log_linear_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
    elif interpolation_method == 'linear':
        print('Starting the interpolation using the "Linear" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            linear_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = linear_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)

    elif interpolation_method == 'log-spline':
        print('Starting the interpolation using the "Log-Spline" method...') if verbose == 1 else None
        interpolated_aep = xarray.apply_ufunc(
            log_spline_interpolate,
            input_flood_wse_backfilled,
            bc_aep,
            elev,
            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes = ['float32'],
            )
        
        if sample_point is not None:
            sample_out = log_spline_interpolate(input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      bc_aep.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      elev.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)

    print('Finished the interpolation') if verbose == 1 else None
    interpolated_aep = interpolated_aep.load()
    interpolated_aep = interpolated_aep.where(input_flood_wse.isel(returnPeriod=n_inputs-1) >= -1000,np.nan)
    
    if extrapolate == 0:
        interpolated_aep = interpolated_aep.where((interpolated_aep <= np.max(raster_exc_probs)) | np.isnan(interpolated_aep),np.max(raster_exc_probs))
        interpolated_aep = interpolated_aep.where((interpolated_aep >= np.min(raster_exc_probs)) | np.isnan(interpolated_aep),np.min(raster_exc_probs))
        
        print(f'Converted AEP values greater than {np.max(raster_exc_probs)} to {np.max(raster_exc_probs)} and values less than {np.min(raster_exc_probs)} to NaN') if verbose == 1 else None
        
    if aep_raster_generate == 1:
        generate_tif(elevation_raster,interpolated_aep,output_data_dir + f'/aep_{interpolation_method}_{areaName}.tif')
        
    if aep_map_generate == 1:
        print('Saving the AEP plot to ...') if verbose == 1 else None
        plt.figure(figsize=(20,15))
        plt.pcolor(interpolated_aep.x,interpolated_aep.y,interpolated_aep,vmin=np.min(raster_exc_probs),vmax=np.max(raster_exc_probs))
        plt.colorbar()
        plt.title('AEP')
        if sample_point is not None:
            plt.scatter(sample_point[0],sample_point[1],color='red')
        plt.savefig(output_plot_dir + f'/aep_map_{interpolation_method}_{areaName}.png')
        plt.close()
        
    if flood_zone_return_period is not None:
        n_output_rp = len(flood_zone_return_period)
        print(f'Generating output for {n_output_rp} new return periods:{flood_zone_return_period}') if verbose == 1 else None
        print('Checking the validity of the return periods...') if verbose == 1 else None
        
        if (np.min(flood_zone_return_period) < np.min(input_return_periods)) and (extrapolate == 0):
            print('Error: extrapolate is set to zero but a return period outside of the input return period range is requested') if verbose == 1 else None
        elif (np.max(flood_zone_return_period) > np.max(input_return_periods)) and (extrapolate == 0):
            print('Error: extrapolate is set to zero but a return period outside of the input return period range is requested') if verbose == 1 else None

        flood_zone_return_period = np.array(flood_zone_return_period)
        output_exc_prob = 1/flood_zone_return_period
        
        for j in range(0,len(flood_zone_return_period)):
            output_rp_zone = interpolated_aep.where(interpolated_aep >= output_exc_prob[j], np.nan)
            
            if flood_zone_raster_generate == 1:
                generate_tif(elevation_raster,output_rp_zone,output_data_dir + f'/estimated_aep_data_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.tif')
                
            dataC = np.zeros(elev.shape,float)
            dataC[:] = flood_zone_return_period[j]
            return_period_da = elev.copy(deep=True,data=dataC)

            # Start the interpolation
            if interpolation_method == 'spline':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Spline method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            spline_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = spline_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None
                    
            elif interpolation_method == 'log-linear':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Log-Linear method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            log_linear_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = log_linear_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

            elif interpolation_method == 'linear':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Linear method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            linear_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = linear_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

            elif interpolation_method == 'log-spline':
                print(f'Starting the interpolation of water surface elevation for return period {flood_zone_return_period[j]} using Log-Spline method...') if verbose == 1 else None
                interpolated_wse = xarray.apply_ufunc(
                            log_spline_interpolate,
                            bc_rp,
                            input_flood_wse_backfilled,
                            return_period_da,
                            input_core_dims=[['returnPeriod'],['returnPeriod'],[]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes = ['float32'],
                            )
                
                if sample_point is not None:
                    sample_out = log_spline_interpolate(bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values,
                                      plot=True)
                    print(f"Received coordinates of a sample point. Return periods: {bc_rp.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n water surface elevation: {input_flood_wse_backfilled.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output return period requested: {return_period_da.sel(x=sample_point[0],y=sample_point[1],method='nearest').values}\n output water surface elevation: {sample_out} ") if verbose == 1 else None

    
            interpolated_wse = interpolated_wse.where(interpolated_aep >= output_exc_prob[j], np.nan)
            generate_tif(elevation_raster,interpolated_wse.values,output_data_dir + f'/estimated_wse_data_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.tif')
            print(f'Finished the interpolation of water surface elevation for the requested new return period: {flood_zone_return_period[j]}') if verbose == 1 else None

            
            if flood_zone_map_generate == 1:
                plt.figure(figsize=(20,15))
                plt.pcolor(output_rp_zone.x,output_rp_zone.y,output_rp_zone,vmin=np.min(output_rp_zone),vmax=np.max(output_rp_zone))
                plt.colorbar()
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood AEP')
                plt.savefig(output_plot_dir + f'/estimated_aep_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
                
                output_rp_aep_inout = output_rp_zone.where(np.isnan(output_rp_zone),1)
                plt.figure(figsize=(20,15))
                plt.pcolor(output_rp_aep_inout.x,output_rp_aep_inout.y,output_rp_aep_inout,vmin=0,vmax=1/100)
                plt.colorbar()
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood in-out zone')
                plt.savefig(output_plot_dir + f'/estimated_aep_zonal_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
                
                plt.figure(figsize=(20,15))
                plt.pcolor(interpolated_wse.x,interpolated_wse.y,interpolated_wse)
                if sample_point is not None:
                    plt.scatter(sample_point[0],sample_point[1],color='red')
                plt.colorbar()
                plt.title(f'Estimated {flood_zone_return_period[j]}-year flood WSE')
                plt.savefig(output_plot_dir + f'/estimated_wse_map_rp_{flood_zone_return_period[j]}_{interpolation_method}_{areaName}.png')
                plt.close()
    return
    
def backfill_nan_with_nearest(darr: xarray.DataArray):
    """
    Given array with member_id/time/lat/lon dimensions, will go through all times in the array and 
    backfill NaN values with nearest non-NaN neighbor
    Assumes static mask in time
    """
    # Part 1...build the NaN Mask
    nmem = len(darr.returnPeriod)
    print(nmem)
    ny = len(darr.y)
    nx = len(darr.x)
    print(f'Received an array with shape {nx,ny,nmem} for backfilling')
    # Build nan mask
    #thisarr = np.ma.masked_where(np.isnan(varray[t,z].values), varray[t,z].values)
    # Now loop through time and update
    outarr = darr.values
    for m in range(nmem):
        print(m)
        raw = np.ma.masked_where(np.isnan(darr.isel(returnPeriod=m)),darr.isel(returnPeriod=m))
        # Get array of indices
        y,x=np.mgrid[0:ny,0:nx]
        xygood = np.array((y[~raw.mask],x[~raw.mask])).T
        xybad = np.array((y[raw.mask],x[raw.mask])).T
        print(xygood.shape,xybad.shape)
        the_nearest = KDTree(xygood).query(xybad)[1]
        thisarr = outarr[m,:,:]
        thisarr[raw.mask] = thisarr[~raw.mask][the_nearest]
        outarr[m,:,:] = thisarr
    # Re-wrap in xarray trappings
    outarr = xarray.DataArray(outarr, dims=darr.dims, coords=darr.coords, attrs=darr.attrs)
    return outarr

