import numpy as np
import netCDF4 as nc
import xarray as xr

ds       =xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_HWSD_soil_layer.nc", decode_times=False)
ds_cluster=xr.open_dataset("/data1/zxy/SOC_data_calibration/global_10_clusters_version2_cable.nc", decode_times=False)

ds = ds.where(ds != -99999, 0)

SoilMoist_new = (ds['SoilMoist'] * ds['PFTfrac']).sum(dim='PFT')
SoilTemp_new=(ds['SoilTemp'] * ds['PFTfrac']).sum(dim='PFT')
B_litter_fal_new = (ds['Belowground_litter_fall'] * ds['PFTfrac']).sum(dim='PFT')
Leaf_fall_new = (ds['Leaf_fall'] * ds['PFTfrac']).sum(dim='PFT')
N_litterfall_new = (ds['non_leaf_aboveground_litterfall'] * ds['PFTfrac']).sum(dim='PFT')





ds_new= xr.Dataset({
    'SoilMoist': SoilMoist_new,
    'SoilTemp': SoilTemp_new,    
    'Belowground_litter_fall': B_litter_fal_new,
    'Leaf_fall': Leaf_fall_new,
    'non_leaf_aboveground_litterfall': N_litterfall_new,
    'water_potential': ds['water_potential']
})

ds_new['soil_cluster'] = xr.DataArray(
    ds_cluster['Band1'].values,
    dims=['lat', 'lon'],
    coords={'lat': ds_new.lat, 'lon': ds_new.lon}
)


ds_new['soil_cluster'].attrs['long_name'] = "global_10_clusters_version2"
ds_new['SoilMoist'].attrs['long_name'] = "Soil moisture for HWSD layers"
ds_new['SoilMoist'].attrs['units'] = "m3/m3"

ds_new['SoilTemp'].attrs['long_name'] = "Soil temperature for HWSD layers"
ds_new['SoilTemp'].attrs['units'] = "degree Celsius(°C)"

ds_new['water_potential'].attrs['long_name'] = "water potential for HWSD layers"
ds_new['water_potential'].attrs['units'] = "kPa"
if 'description' in ds_new['water_potential'].attrs:
    del ds_new['water_potential'].attrs['description']

ds_new = ds_new.where(ds_new != 0, -99999)
ds_new = ds_new.fillna(-99999)

ds_new = ds_new.drop_vars('time', errors='ignore')

ds_new.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_avePFT.nc")