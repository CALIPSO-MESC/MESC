import numpy as np
import netCDF4 as nc
import xarray as xr



ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/data_extraction_cable_model_calibration_cluster.nc", decode_times=False)
dataset_invariant=xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

HWSD_SOC = dataset_invariant['HWSD_SOC']   
HWSD_bulk_density = dataset_invariant['HWSD_bulk_density']   

HWSD_SOC_set=[]
HWSD_bulk_density_set=[]
for lat, lon in zip(ds['lat'].values, ds['lon'].values):
    print(lat,lon)
    soc_value = HWSD_SOC.sel(lat=lat, lon=lon)
    bd_value = HWSD_bulk_density.sel(lat=lat, lon=lon)
    print(len(soc_value),len(bd_value))
    print(soc_value.values)
    print(bd_value.values)

    HWSD_SOC_set.append(soc_value)
    HWSD_bulk_density_set.append(bd_value)

HWSD_SOC_set = np.array(HWSD_SOC_set)
HWSD_bulk_density_set = np.array(HWSD_bulk_density_set)
print(HWSD_SOC_set.shape)

HWSD_SOC_set = HWSD_SOC_set.T
HWSD_bulk_density_set = HWSD_bulk_density_set.T
print(HWSD_SOC_set.shape)

ds['HWSD_SOC'].values=HWSD_SOC_set

ds['HWSD_bulk_density'] = ds['HWSD_SOC'].copy()
ds['HWSD_bulk_density'].values=HWSD_bulk_density_set
ds['HWSD_bulk_density'].attrs['long_name'] = "HWSD2 Soil bulk density"
ds['HWSD_bulk_density'].attrs['units'] = "kg/m3"

ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/data_extraction_cable_model_calibration_cluster_new.nc")
