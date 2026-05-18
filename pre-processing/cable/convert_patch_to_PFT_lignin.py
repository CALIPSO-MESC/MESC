import numpy as np
import netCDF4 as nc


dataset=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/lignin_C_cable.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/cable_out_hwsd90_cn_1990_extracted_var.nc")

# SoilTemp(time, soil, patch, y, x) 
lignin_leaf_grid=dataset['lignin_leaf_grid'][:]
lignin_CWD_grid=dataset['lignin_CWD_grid'][:]
lignin_froot_grid=dataset['lignin_froot_grid'][:]

patchfrac=dataset_static['patchfrac'][:]
iveg=dataset_static['iveg'][:]

lignin_leaf_new  = np.full((17,112,192), np.nan, dtype=float)
lignin_CWD_new  = np.full((17,112,192), np.nan, dtype=float)
lignin_froot_new  = np.full((17,112,192), np.nan, dtype=float)

for PFT in range(1,18):
    print(PFT)    
    lignin_leaf_new[PFT-1,:,:] = np.nanmean(np.where(iveg == PFT, lignin_leaf_grid, np.nan),axis=0)
    lignin_CWD_new[PFT-1,:,:] = np.nanmean(np.where(iveg == PFT, lignin_CWD_grid, np.nan),axis=0)
    lignin_froot_new[PFT-1,:,:] = np.nanmean(np.where(iveg == PFT, lignin_froot_grid, np.nan),axis=0)    

new = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/lignin_PFT.nc", 'w')

new.createDimension("lat", 112 )
new.createDimension("lon", 192 )
new.createDimension("PFT", 17 )

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=dataset_static['y'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=dataset_static['x'][:]

new.createVariable("PFT", "float64", (u'PFT',))
new.variables['PFT'].__setattr__("long_name","plant function type")
new.variables['PFT'][:]=range(1,18)


new.createVariable("lignin_leaf", "float64", (u'PFT',u'lat',u'lon',))
new.variables['lignin_leaf'].__setattr__("long_name","Lignin:C ratio leaf litter fall")
new.variables['lignin_leaf'][:]=lignin_leaf_new[:]

new.createVariable("lignin_CWD", "float64", (u'PFT',u'lat',u'lon',))
new.variables['lignin_CWD'].__setattr__("long_name","Lignin:C ratio non-leaf litter fall")
new.variables['lignin_CWD'][:]=lignin_CWD_new[:]

new.createVariable("lignin_froot", "float64", (u'PFT',u'lat',u'lon',))
new.variables['lignin_froot'].__setattr__("long_name","Lignin:C ratio of belowground litter fall")
new.variables['lignin_froot'][:]=lignin_froot_new[:]


