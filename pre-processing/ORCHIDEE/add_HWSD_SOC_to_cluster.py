import netCDF4 as nc
import pandas as pd
import numpy as np
import xarray as xr
import sys
import os


#===================================write nc file========================================#

dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/data_extraction_ORCHIDEE_model_calibration.nc", decode_times=False)
ds= xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/data_extraction_ORCHIDEE_model_calibration_v2_cluster.nc", decode_times=False)

ds['HWSD_SOC'][:]=dataset['HWSD_SOC'].values


ds['HWSD_bulk_density'] = ds['HWSD_SOC'].copy()
ds['HWSD_bulk_density'][:]=dataset['HWSD_bulk_density'].values
ds['HWSD_bulk_density'].attrs['long_name'] = "HWSD2 Soil bulk density"
ds['HWSD_bulk_density'].attrs['units'] = "kg/m3"


ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/data_extraction_ORCHIDEE_model_calibration_v2_cluster_new.nc")





