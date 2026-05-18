#convert ORCHIDEE soil layer to HWSD soil layer

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
import sys

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress_new.nc", decode_times=False)
SoilTemp_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/soiltemp_daily/soiltemp1990_2009_new.nc", decode_times=False)
SoilMoist_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/soilmoist_daily/ORCHIDEE_soilmoist_ave_1990_2020_time_varing_new.nc", decode_times=False)
water_potential_data=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/water_potential_USDA_Texture_ORCHIDEE_new.nc", decode_times=False)

ds['Leaf_fall'] = xr.where(ds['Leaf_fall'] > 10000000, -99999, ds['Leaf_fall'])
ds['Belowground_litter_fall'] = xr.where(ds['Belowground_litter_fall'] > 10000000, -99999, ds['Belowground_litter_fall'])
ds['non_leaf_aboveground_litterfall'] = xr.where(ds['non_leaf_aboveground_litterfall'] > 10000000, -99999, ds['non_leaf_aboveground_litterfall'])

#Soiltemp(time, solth, lat, lon)  soilmoist(time, veget, solay, lat, lon)
SoilTemp = SoilTemp_data['Soiltemp'].values-273.15
SoilMoist=SoilMoist_data['soilmoist'].values
water_potential=water_potential_data['water_potential'].values
print("load data")
sys.stdout.flush()

ntime,  npft, nsoil_orig, nlat, nlon = SoilMoist.shape
nsoil_new = 7

SoilTemp_HWSD = xr.DataArray(
    np.zeros((ntime, nsoil_new, nlat, nlon)),
    dims=('time',  'soil_hwsd',  'lat', 'lon'),
    coords={'time': ds['time'], 'soil_hwsd': np.arange(nsoil_new),'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'Soil temperature for HWSD layers',
        'units': 'degree Celsius(°C)',
        'description': 'Soil temperature with HWSD soil layers'
    }
)

SoilMoist_HWSD = xr.DataArray(
    np.zeros((ntime, npft, nsoil_new, nlat, nlon)),
    dims=('time', 'veget', 'soil_hwsd',  'lat', 'lon'),
    coords={'time': ds['time'],'veget': ds['veget'], 'soil_hwsd': np.arange(nsoil_new),'lat': ds['lat'], 'lon': ds['lon']},
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

print("star cal")
sys.stdout.flush()

SoilTemp_HWSD[:,0,:,:]=(0.097752 *SoilTemp[:,0,:,:] + 0.293255 *SoilTemp[:,1,:,:]
                       +0.586510 *SoilTemp[:,2,:,:] + 1.173021 *SoilTemp[:,3,:,:]
                       +2.346041 *SoilTemp[:,4,:,:] + 4.692082 *SoilTemp[:,5,:,:]
                       +9.384164 *SoilTemp[:,6,:,:] + 1.427175 *SoilTemp[:,7,:,:])/20
SoilTemp_HWSD[:,1,:,:]=(17.341153*SoilTemp[:,7,:,:] + 2.658847 *SoilTemp[:,8,:,:])/20
SoilTemp_HWSD[:,2,:,:]=SoilTemp[:,8,:,:]
SoilTemp_HWSD[:,3,:,:]=(14.877810*SoilTemp[:,8,:,:] + 5.122190 *SoilTemp[:,9,:,:])/20
SoilTemp_HWSD[:,4,:,:]=SoilTemp[:,9,:,:]
SoilTemp_HWSD[:,5,:,:]=(49.951124*SoilTemp[:,9,:,:] + 0.048876 *SoilTemp[:,10,:,:])/50
SoilTemp_HWSD[:,6,:,:]=SoilTemp[:,10,:,:]
print("SoilTemp_HWSD")
sys.stdout.flush()

SoilMoist_HWSD[:,:,0,:,:]=(0.097752 *SoilMoist[:,:,0,:,:] + 0.293255 *SoilMoist[:,:,1,:,:]
                          +0.586510 *SoilMoist[:,:,2,:,:] + 1.173021 *SoilMoist[:,:,3,:,:]
                          +2.346041 *SoilMoist[:,:,4,:,:] + 4.692082 *SoilMoist[:,:,5,:,:]
                          +9.384164 *SoilMoist[:,:,6,:,:] + 1.427175 *SoilMoist[:,:,7,:,:])/20
SoilMoist_HWSD[:,:,1,:,:]=(17.341153*SoilMoist[:,:,7,:,:] + 2.658847 *SoilMoist[:,:,8,:,:])/20
SoilMoist_HWSD[:,:,2,:,:]=SoilMoist[:,:,8,:,:]
SoilMoist_HWSD[:,:,3,:,:]=(14.877810*SoilMoist[:,:,8,:,:] + 5.122190 *SoilMoist[:,:,9,:,:])/20
SoilMoist_HWSD[:,:,4,:,:]=SoilMoist[:,:,9,:,:]
SoilMoist_HWSD[:,:,5,:,:]=(49.951124*SoilMoist[:,:,9,:,:] + 0.048876 *SoilMoist[:,:,10,:,:])/50
SoilMoist_HWSD[:,:,6,:,:]=SoilMoist[:,:,10,:,:]
print("SoilMoist_HWSD")
sys.stdout.flush()

water_potential_HWSD[:,0,:,:]=(0.097752 *water_potential[:,0,:,:] + 0.293255 *water_potential[:,1,:,:]
                              +0.586510 *water_potential[:,2,:,:] + 1.173021 *water_potential[:,3,:,:]
                              +2.346041 *water_potential[:,4,:,:] + 4.692082 *water_potential[:,5,:,:]
                              +9.384164 *water_potential[:,6,:,:] + 1.427175 *water_potential[:,7,:,:])/20
water_potential_HWSD[:,1,:,:]=(17.341153*water_potential[:,7,:,:] + 2.658847 *water_potential[:,8,:,:])/20
water_potential_HWSD[:,2,:,:]=water_potential[:,8,:,:]
water_potential_HWSD[:,3,:,:]=(14.877810*water_potential[:,8,:,:] + 5.122190 *water_potential[:,9,:,:])/20
water_potential_HWSD[:,4,:,:]=water_potential[:,9,:,:]
water_potential_HWSD[:,5,:,:]=(49.951124*water_potential[:,9,:,:] + 0.048876 *water_potential[:,10,:,:])/50
water_potential_HWSD[:,6,:,:]=water_potential[:,10,:,:]

print("write nc file")
sys.stdout.flush()


# 替换SoilMoist和water_potential
ds['SoilTemp'].values = np.where(ds['SoilTemp'].values != -99999, SoilTemp_HWSD.values, ds['SoilTemp'].values)
ds['SoilMoist'].values = np.where(ds['SoilMoist'].values != -99999, SoilMoist_HWSD.values, ds['SoilMoist'].values)
ds['water_potential'].values = np.where(ds['water_potential'].values != -99999, water_potential_HWSD.values, ds['water_potential'].values)


# 保存到新文件
ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress_new_new.nc")


