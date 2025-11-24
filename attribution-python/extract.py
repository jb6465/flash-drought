'''
############################################################################################
## BUREAU OF METEOROLOGY
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT AND CHANGES IN ARIDITY HAZARD TEAM
##
## DATE:             Nov-2025
## SCRIPT:           extract.py
## AUTHOR:           jessica.bhardwaj@bom.gov.au
##
## DESCRIPTION:      Script to extract reanalysis and projections datasets.
##
############################################################################################
'''

# directories
barraR2_dir = '/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/'
era5_dir = '/g/data/rt52/era5/'
merra2_M2T1NXSLV_dir = '/g/data/rr7/MERRA2/raw/M2T1NXSLV.5.12.4'
merra2_M2T1NXRAD_dir = '/g/data/rr7/MERRA2/raw/M2T1NXRAD.5.12.4/'

domain_dict = { #can add other domains for future work
    'australia':{'lat_min':-44, 'lat_max':-10, 'lon_min':112, 'lon_max':154},
    'barra_domain':{'lat_min':-57.97, 'lat_max':12.98, 'lon_min':88.48, 'lon_max':207.39},
}

#functions
def barra_daily_extract(target_var, extracted_data_save_dir, nation_domain, year):
    """
    Inputs:
    - target_var - string of BARRA-R2 target variable key
    - extracted_data_save_dir - string of directory to save extracted data in
    - nation_domain - string to specify target nation domain 
    - year - string of target extraction year
    Returns:
    Daily BARRA variable files downloaded and saved in specified dir.
    """
    
    import os
    import glob
    import xarray as xr
    
    if not os.path.isfile(f"{extracted_data_save_dir}/{target_var}/{nation_domain}_BARRAR2_{year}_{target_var}_day.nc"):
        os.makedirs(os.path.abspath(f"{extracted_data_save_dir}/{target_var}/"), exist_ok=True)
        barra_files = sorted(glob.glob(f"{barraR2_dir}/{target_var}/latest/*{year}*.nc"))
        datasets = []
        for file in barra_files:
            ds = xr.open_dataset(file, engine='netcdf4')
            datasets.append(ds[target_var].sel(lat=slice(domain_dict[nation_domain]['lat_min'], domain_dict[nation_domain]['lat_max']), lon=slice(domain_dict[nation_domain]['lon_min'], domain_dict[nation_domain]['lon_max'])))
            del ds
        ds_cube = xr.concat(datasets, dim='time')
        ds_cube.to_netcdf(f"{extracted_data_save_dir}/{target_var}/{nation_domain}_BARRAR2_{year}_{target_var}_day.nc", encoding={f"{target_var}": {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
        ds_cube.close()
    return

def era5_daily_extract(target_var, extracted_data_save_dir, nation_domain, year):
    """
    Inputs:
    - target_var - string of ERA5 target variable key
    - extracted_data_save_dir - string of directory to save extracted data in
    - nation_domain - string to specify target nation domain 
    - year - string of target extraction year
    Returns:
    Hourly ERA5 variable files downloaded, aggregated to daily and saved in specified dir.
    """
    import os
    import glob
    import xarray as xr
    
    def preprocess_era5(ds):
            ds = ds.rename({'latitude': 'lat'}).rename({'longitude': 'lon'})
            return ds.sel(lat=slice(domain_dict[nation_domain]['lat_max'], domain_dict[nation_domain]['lat_min']), lon=slice(domain_dict[nation_domain]['lon_min'], domain_dict[nation_domain]['lon_max']))
        
    if not os.path.isfile(f"{extracted_data_save_dir}/{target_var}/{nation_domain}_ERA5_{year}_{target_var}_day.nc"):
        os.makedirs(f"{extracted_data_save_dir}/{target_var}/", exist_ok=True)   
        if target_var == '10w':
            u10_files = sorted(glob.glob(f"{era5_dir}single-levels/reanalysis/10u/{year}/*.nc"))
            v10_files = sorted(glob.glob(f"{era5_dir}single-levels/reanalysis/10v/{year}/*.nc"))
            u10_cube = xr.open_mfdataset(u10_files,combine='by_coords',parallel=True, engine='netcdf4', preprocess=preprocess_era5).chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
            v10_cube = xr.open_mfdataset(v10_files,combine='by_coords',parallel=True, engine='netcdf4', preprocess=preprocess_era5).chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
        
            w10 = (v10_cube.v10**2 + u10_cube.u10**2) ** 0.5
            w10 = (w10.sortby("time")).resample(time='D').mean()
            w10 = w10.chunk({'time':720, 'lat':'auto', 'lon':'auto'}).compute()
            w10_cube = xr.Dataset({'w10': w10})
            w10_cube.to_netcdf(f"{extracted_data_save_dir}/{target_var}/{nation_domain}_ERA5_{year}_{target_var}_day.nc", encoding={f"{list(w10_cube.data_vars)[0]}": {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            w10_cube.close()
        else:
            era5_files = sorted(glob.glob(f"{era5_dir}single-levels/reanalysis/{target_var}/{year}/*.nc")) 
            ds_cube = xr.open_mfdataset(era5_files,combine='by_coords',parallel=True, engine='netcdf4', preprocess=preprocess_era5).chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
            ds_cube = (ds_cube.sortby("time")).resample(time='D').mean()
            ds_cube = ds_cube.chunk({'time':720, 'lat':'auto', 'lon':'auto'}).compute()
            ds_cube.to_netcdf(f"{extracted_data_save_dir}/{target_var}/{nation_domain}_ERA5_{year}_{target_var}_day.nc", encoding={f"{list(ds_cube.data_vars)[0]}": {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            ds_cube.close()
    return

def merra2_daily_extract(target_var, extracted_data_save_dir, nation_domain, year):
    """
    Inputs:
    - target_var - string of MERRA2 target variable keys
    - extracted_data_save_dir - string of directory to save extracted data in
    - nation_domain - string to specify target nation domain
    - year - string of target extraction year
    Returns:
    Hourly MERRA2 variable files downloaded, aggregated to daily and saved in specified dir
    """
    import os
    import glob
    import xarray as xr
    import gc
    import warnings
    import logging
    warnings.filterwarnings('ignore') 
    logging.getLogger("distributed").setLevel(logging.ERROR)
    logging.getLogger('flox').setLevel(logging.WARNING)
    
    #preprocess functions to save memory and time
    def preprocess_U2M_V2M(ds):
        logging.getLogger('flox').setLevel(logging.WARNING)
        return ds[['V2M', 'U2M']].sel(lat=slice(domain_dict[nation_domain]['lat_min'], domain_dict[nation_domain]['lat_max']), lon=slice(domain_dict[nation_domain]['lon_min'], domain_dict[nation_domain]['lon_max']))

    def preprocess_T2M_QV2M_PS(ds):
        logging.getLogger('flox').setLevel(logging.WARNING)
        return ds[['T2M', 'QV2M', 'PS']].resample(time='1D').mean().sel(lat=slice(domain_dict[nation_domain]['lat_min'], domain_dict[nation_domain]['lat_max']), lon=slice(domain_dict[nation_domain]['lon_min'], domain_dict[nation_domain]['lon_max']))

    def preprocess_SWdn(ds): 
        logging.getLogger('flox').setLevel(logging.WARNING)
        return ds['SWGDN'].resample(time='1D').mean().sel(lat=slice(domain_dict[nation_domain]['lat_min'], domain_dict[nation_domain]['lat_max']), lon=slice(domain_dict[nation_domain]['lon_min'], domain_dict[nation_domain]['lon_max']))
    
    if target_var == ['U2M', 'V2M']:
        if not os.path.isfile(f"{extracted_data_save_dir}/W2M/{nation_domain}_MERRA2_{year}_W2M_day.nc"):
            files = sorted(glob.glob(f"{merra2_M2T1NXSLV_dir}/{year}/*/*.nc4"))
            U2M_V2M_cube = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True, preprocess=preprocess_U2M_V2M, engine='netcdf4').chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
            U2M_V2M_cube.to_netcdf(f'{extracted_data_save_dir}/U2M_V2M/{nation_domain}_{year}_MERRA2_hly_U2M_V2M.nc', encoding={'U2M': {'zlib': True, 'complevel': 5, 'dtype':'float32'}, 'V2M': {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            del U2M_V2M_cube; gc.collect()
    
            wind = xr.open_mfdataset([f'{extracted_data_save_dir}/U2M_V2M/{file}' for file in os.listdir(f'{extracted_data_save_dir}/U2M_V2M/') if str(year) in file], combine='nested', concat_dim='time', parallel=True, engine='netcdf4').chunk({'time': -1, 'lat': 'auto', 'lon': 'auto'})
            W2M = ((((wind['U2M']**2)+(wind['V2M']**2))**0.5).rename('W2M')).resample(time='1D').mean()
            W2M.to_netcdf(f'{extracted_data_save_dir}/W2M/{nation_domain}_MERRA2_{year}_W2M_day.nc', encoding={'W2M': {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            del W2M; gc.collect()

    if target_var == ['T2M', 'QV2M', 'PS']:
        if not os.path.isfile(f"{extracted_data_save_dir}/T2M_QV2M/{nation_domain}_MERRA2_{year}_T2M_QV2M_day.nc"):
            files = sorted(glob.glob(f"{merra2_M2T1NXSLV_dir}/{year}/*/*.nc4"))
            T2M_QV2M_cube = xr.open_mfdataset(files, combine='nested', concat_dim='time', parallel=True, preprocess=preprocess_T2M_QV2M_PS, engine='netcdf4').chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
            T2M_QV2M_cube.to_netcdf(f'{extracted_data_save_dir}/T2M_QV2M_PS/{nation_domain}_MERRA2_{year}_T2M_QV2M_PS_day.nc', encoding={'T2M': {'zlib': True, 'complevel': 5, 'dtype':'float32'}, 'QV2M': {'zlib': True, 'complevel': 5, 'dtype':'float32'}, 'PS': {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            del T2M_QV2M_cube; gc.collect()

    if target_var == ['SWGDN']:
        if not os.path.isfile(f"{extracted_data_save_dir}/SWGDN/{nation_domain}_MERRA2_{year}_SWGDN_day.nc"):
            files = sorted(glob.glob(f"{merra2_M2T1NXRAD_dir}/{year}/*/*.nc4"))
            SWGDN_cube = xr.open_mfdataset(files,combine='nested', concat_dim='time',parallel=True, preprocess=preprocess_SWdn, engine='netcdf4').chunk({'time':-1, 'lat':'auto', 'lon':'auto'})
            SWGDN_cube.to_netcdf(f'{extracted_data_save_dir}/SWGDN/{nation_domain}_MERRA2_{year}_SWGDN_day.nc', encoding={'SWGDN': {'zlib': True, 'complevel': 5, 'dtype':'float32'}})
            del SWGDN_cube; gc.collect()
    
    return