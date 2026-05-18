import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
import sys

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress.nc", decode_times=False)
litter_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/litter_daily/ORCHIDEE_litter_ave_1990_2009_time_varing.nc", decode_times=False)

#Soiltemp(time, solth, lat, lon)  soilmoist(time, veget, solay, lat, lon)
Leaf_fall = litter_data['Leaf_fall'].values
non_leaf_aboveground_litterfall=litter_data['non_leaf_aboveground_litterfall'].values
Belowground_litter_fall=litter_data['Belowground_litter_fall'].values
print("load data")
sys.stdout.flush()


ds['Leaf_fall'].values  = Leaf_fall
ds['non_leaf_aboveground_litterfall'].values  = non_leaf_aboveground_litterfall
ds['Belowground_litter_fall'].values  = Belowground_litter_fall


# 保存到新文件
ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress_new.nc")


