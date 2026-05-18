import numpy as np
import netCDF4 as nc
import xarray as xr

ds_static=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)
ds       =xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress.nc", decode_times=False)
ds_cluster=xr.open_dataset("/data1/zxy/SOC_data_calibration/global_10_clusters_version2.nc", decode_times=False)

ds_static = ds_static.squeeze('time', drop=True)

ds = ds.where(ds != -99999, 0)
ds_static = ds_static.where(ds_static != -99999, 0)

SoilMoist_new = (ds['SoilMoist'] * ds_static['maxvegetfrac']).sum(dim='veget')
SoilTemp_new = (ds['SoilTemp'] * ds_static['maxvegetfrac']).sum(dim='veget')
B_litter_fal_new = (ds['Belowground_litter_fall'] * ds_static['maxvegetfrac']).sum(dim='veget')
Leaf_fall_new = (ds['Leaf_fall'] * ds_static['maxvegetfrac']).sum(dim='veget')
N_litterfall_new = (ds['non_leaf_aboveground_litterfall'] * ds_static['maxvegetfrac']).sum(dim='veget')

ds_new= xr.Dataset({
    'SoilMoist': SoilMoist_new,
    'SoilTemp': SoilTemp_new,
    'Belowground_litter_fall': B_litter_fal_new,
    'Leaf_fall': Leaf_fall_new,
    'non_leaf_aboveground_litterfall': N_litterfall_new,
    'water_potential': ds['water_potential']
})

ds_new['soil_cluster'] = xr.DataArray(
    ds_cluster['Band1'].values[::-1, :],
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

ds_new.to_netcdf('/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_PTFave.nc')