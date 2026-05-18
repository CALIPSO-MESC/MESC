import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
import sys

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)
pH_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/pH_0-20cm_ave_soilgrid_0.5d.nc", decode_times=False)
clay_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/clay_0-20cm_ave_soilgrid_0.5d.nc", decode_times=False)
silt_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/silt_0-20cm_ave_soilgrid_0.5d.nc", decode_times=False)

SoilSuborder_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_0.5d.nc", decode_times=False)


ORCHIDEE_Soil_texture_class=np.array(ds['Soil_texture'])

clay=np.array([0.03,	0.06,	0.11,	0.19,	0.1,	0.2,	0.27,	0.33,	0.33,	0.41,	0.46,	0.55])
silt=np.array([0.04,	0.13,	0.26,	0.64,	0.84,	0.4,	0.19,	0.59,	0.37,	0.11,	0.48,	0.3])


clay_ORCHIDEE=ORCHIDEE_Soil_texture_class.copy()
silt_ORCHIDEE=ORCHIDEE_Soil_texture_class.copy()


mask = ORCHIDEE_Soil_texture_class != -99999

clay_ORCHIDEE[mask] = clay[ORCHIDEE_Soil_texture_class[mask].astype(int) - 1]
silt_ORCHIDEE[mask] = silt[ORCHIDEE_Soil_texture_class[mask].astype(int) - 1]



ds['clay_ORCHIDEE'] = ds['Ald'].copy()
ds['silt_ORCHIDEE'] = ds['Ald'].copy()


ds['clay_ORCHIDEE'].attrs = {}
ds['silt_ORCHIDEE'].attrs = {}



ds['clay_ORCHIDEE'].values  = clay_ORCHIDEE
ds['silt_ORCHIDEE'].values  = silt_ORCHIDEE



ds['pH_soilgrid'] = ds['Ald'].copy()
ds['clay_soilgrid'] = ds['Ald'].copy()
ds['silt_soilgrid'] = ds['Ald'].copy()
ds['USDA_SoilSuborder'] = ds['Ald'].copy()

ds['pH_soilgrid'].attrs = {}
ds['clay_soilgrid'].attrs = {}
ds['silt_soilgrid'].attrs = {}
ds['USDA_SoilSuborder'].attrs = {}

ds['pH_soilgrid'].values  = pH_data['pH0_20cm'][::-1, :].where(lambda x: x != 0, -99999)
ds['clay_soilgrid'].values  = clay_data['clay0_20cm'][::-1, :].where(lambda x: x != 0, -99999)
ds['silt_soilgrid'].values  = silt_data['silt0_20cm'][::-1, :].where(lambda x: x != 0, -99999)
ds['USDA_SoilSuborder'].values  = SoilSuborder_data['USDA_soil_map'][::-1, :].fillna(-99999)

ds['pH_soilgrid'].attrs['depth'] = "0-20cm"
ds['clay_soilgrid'].attrs['depth'] = "0-20cm"
ds['silt_soilgrid'].attrs['depth'] = "0-20cm"

ds['pH_soilgrid'].attrs['units'] = "g/100g(%)"
ds['clay_soilgrid'].attrs['units'] = "g/100g(%)"
ds['silt_soilgrid'].attrs['units'] = "g/100g(%)"

ds['USDA_SoilSuborder'].attrs['long_name'] = "USDA Soil texture class"



# 保存到新文件
ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new.nc")