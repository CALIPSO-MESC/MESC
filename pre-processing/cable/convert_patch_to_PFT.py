import numpy as np
import netCDF4 as nc


dataset=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/hwsd90_cn_2990-3009_ave_daily_SoilTemp_SoilMoist.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/cable_out_hwsd90_cn_1990_extracted_var.nc")

# SoilTemp(time, soil, patch, y, x) 
SoilTemp=dataset['SoilTemp'][:]
SoilMoist=dataset['SoilMoist'][:]

patchfrac=dataset_static['patchfrac'][:]
iveg=dataset_static['iveg'][:]

SoilTemp_new  = np.full((365,6,17,112,192), np.nan, dtype=float)
SoilMoist_new  = np.full((365,6,17,112,192), np.nan, dtype=float)
PFTfrac = np.full((17,112,192), np.nan, dtype=float)

for PFT in range(1,18):
    print(PFT)
    PFTfrac[PFT-1,:,:] = np.nanmean(np.where(iveg == PFT, patchfrac, np.nan),axis=0)
    SoilTemp_new[:,:,PFT-1,:,:] = np.nanmean(np.where(iveg[np.newaxis, np.newaxis, :, :, :] == PFT, SoilTemp, np.nan),axis=2)
    SoilMoist_new[:,:,PFT-1,:,:] = np.nanmean(np.where(iveg[np.newaxis, np.newaxis, :, :, :] == PFT, SoilMoist, np.nan),axis=2)
    

new = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/hwsd90_cn_2990-3009_ave_daily_SoilTemp_SoilMoist_PFT.nc", 'w')

new.createDimension("time", 365 )
new.createDimension("soil", 6 )
new.createDimension("lat", 112 )
new.createDimension("lon", 192 )
new.createDimension("PFT", 17 )


new.createVariable("time", "float64", (u'time',))
new.variables['time'].__setattr__("units","seconds since 1990-01-01 00:00:00")
new.variables['time'].__setattr__("coordinate","GMT")
new.variables['time'].__setattr__("calendar","")
new.variables['time'][:]=dataset['time'][:]

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=dataset_static['y'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=dataset_static['x'][:]

new.createVariable("PFT", "float64", (u'PFT',))
new.variables['PFT'].__setattr__("long_name","plant function type")
new.variables['PFT'][:]=range(1,18)


new.createVariable("SoilMoist", "float64", (u'time',u'soil',u'PFT',u'lat',u'lon',))
new.variables['SoilMoist'].__setattr__("units","m^3/m^3")
new.variables['SoilMoist'].__setattr__("long_name","Average layer soil moisture")
new.variables['SoilMoist'][:]=SoilMoist_new[:]

new.createVariable("SoilTemp", "float64", (u'time',u'soil',u'PFT',u'lat',u'lon',))
new.variables['SoilTemp'].__setattr__("units","K")
new.variables['SoilTemp'].__setattr__("long_name","Average layer soil temperature")
new.variables['SoilTemp'][:]=SoilTemp_new[:]

new.createVariable("PFTfrac", "float64", (u'PFT',u'lat',u'lon',))
new.variables['PFTfrac'].__setattr__("long_name","Fractional cover of PFT")
new.variables['PFTfrac'][:]=PFTfrac[:]

