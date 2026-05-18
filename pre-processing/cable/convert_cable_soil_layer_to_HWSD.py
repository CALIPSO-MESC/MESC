#convert cable soil layer to HWSD soil layer

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables.nc", decode_times=False)

#(time, soil, PFT, lat, lon) 
SoilTemp = ds['SoilTemp'].data-273.15
SoilMoist=ds['SoilMoist'].data
water_potential=ds['water_potential'].data

ntime, nsoil_orig, npft, nlat, nlon = SoilTemp.shape
nsoil_new = 7

SoilTemp_HWSD = xr.DataArray(
    np.zeros((ntime, nsoil_new, npft, nlat, nlon)),
    dims=('time', 'soil_hwsd', 'PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'soil_hwsd': np.arange(nsoil_new),'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'Soil temperature for HWSD layers',
        'units': 'degree Celsius(°C)',
        'description': 'Soil temperature with HWSD soil layers'
    }
)

SoilMoist_HWSD = xr.DataArray(
    np.zeros((ntime, nsoil_new, npft, nlat, nlon)),
    dims=('time', 'soil_hwsd', 'PFT', 'lat', 'lon'),
    coords=SoilTemp_HWSD.coords,
    attrs={
        'long_name': 'Soil moisture for HWSD layers',
        'units': 'm3/m3',
        'description': 'Soil moisture of HWSD soil layers'
    }
)

water_potential_HWSD = xr.DataArray(
    np.zeros((ntime, nsoil_new,  nlat, nlon)),
    dims=('time', 'soil_hwsd',  'lat', 'lon'),
    coords={'time': ds['time'],'soil_hwsd': np.arange(nsoil_new), 'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'water potential for HWSD layers',
        'units': 'kPa',
        'description': 'Soil moisture of HWSD soil layers'
    }
)

SoilTemp[SoilTemp==-99999]=np.nan
SoilMoist[SoilMoist==-99999]=np.nan
water_potential[water_potential==-99999]=np.nan

SoilTemp_HWSD[:,0,:,:,:]=(2.2*SoilTemp[:,0,:,:,:]+5.6*SoilTemp[:,1,:,:,:]+12.2*SoilTemp[:,2,:,:,:])/20
SoilTemp_HWSD[:,1,:,:,:]=(3.2*SoilTemp[:,2,:,:,:]+16.8*SoilTemp[:,3,:,:,:])/20
SoilTemp_HWSD[:,2,:,:,:]=SoilTemp[:,3,:,:,:]
SoilTemp_HWSD[:,3,:,:,:]=(4.1*SoilTemp[:,3,:,:,:]+15.9*SoilTemp[:,4,:,:,:])/20
SoilTemp_HWSD[:,4,:,:,:]=SoilTemp[:,4,:,:,:]
SoilTemp_HWSD[:,5,:,:,:]=SoilTemp[:,4,:,:,:]
SoilTemp_HWSD[:,6,:,:,:]=(22.6*SoilTemp[:,4,:,:,:]+27.4*SoilTemp[:,5,:,:,:])/50

SoilMoist_HWSD[:,0,:,:,:]=(2.2*SoilMoist[:,0,:,:,:]+5.6*SoilMoist[:,1,:,:,:]+12.2*SoilMoist[:,2,:,:,:])/20
SoilMoist_HWSD[:,1,:,:,:]=(3.2*SoilMoist[:,2,:,:,:]+16.8*SoilMoist[:,3,:,:,:])/20
SoilMoist_HWSD[:,2,:,:,:]=SoilMoist[:,3,:,:,:]
SoilMoist_HWSD[:,3,:,:,:]=(4.1*SoilMoist[:,3,:,:,:]+15.9*SoilMoist[:,4,:,:,:])/20
SoilMoist_HWSD[:,4,:,:,:]=SoilMoist[:,4,:,:,:]
SoilMoist_HWSD[:,5,:,:,:]=SoilMoist[:,4,:,:,:]
SoilMoist_HWSD[:,6,:,:,:]=(22.6*SoilMoist[:,4,:,:,:]+27.4*SoilMoist[:,5,:,:,:])/50

water_potential_HWSD[:,0,:,:]=(2.2*water_potential[:,0,:,:]+5.6*water_potential[:,1,:,:]+12.2*water_potential[:,2,:,:])/20
water_potential_HWSD[:,1,:,:]=(3.2*water_potential[:,2,:,:]+16.8*water_potential[:,3,:,:])/20
water_potential_HWSD[:,2,:,:]=water_potential[:,3,:,:]
water_potential_HWSD[:,3,:,:]=(4.1*water_potential[:,3,:,:]+15.9*water_potential[:,4,:,:])/20
water_potential_HWSD[:,4,:,:]=water_potential[:,4,:,:]
water_potential_HWSD[:,5,:,:]=water_potential[:,4,:,:]
water_potential_HWSD[:,6,:,:]=(22.6*water_potential[:,4,:,:]+27.4*water_potential[:,5,:,:])/50



# 替换SoilTemp和SoilMoist变量
ds = ds.drop_vars(['SoilTemp', 'SoilMoist', 'water_potential'])
ds['SoilTemp'] = SoilTemp_HWSD
ds['SoilMoist'] = SoilMoist_HWSD
ds['water_potential'] = water_potential_HWSD

# 添加新的soil变量

ds['soil_hwsd'] = (('soil_hwsd',), np.arange(1,8))
ds['soil_hwsd'].attrs['long_name'] = 'HWSD Soil layer' 
ds['soil_hwsd'].attrs['layer1'] = "0-20cm"
ds['soil_hwsd'].attrs['layer2'] = "20-40cm"
ds['soil_hwsd'].attrs['layer3'] = "40-60cm"
ds['soil_hwsd'].attrs['layer4'] = "60-80cm"
ds['soil_hwsd'].attrs['layer5'] = "80-100cm"
ds['soil_hwsd'].attrs['layer6'] = "100-150cm"
ds['soil_hwsd'].attrs['layer7'] = "150-200cm"

# 保存到新文件
ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_HWSD_soil_layer.nc")


