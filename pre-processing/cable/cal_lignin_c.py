import numpy as np
import pandas as pd
import netCDF4 as nc

dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/cable_out_hwsd90_cn_1990_extracted_var.nc")

lignin_leaf  = np.array([0.25,	0.2,	0.2,	0.2,	0.2,	0.1,	0.1,	0.1,	0.1,	0.1,	0.25,	0.2,	0.2,	0.15,	0.15,	0.25,	0.1])
lignin_CWD   = np.array([0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4,	0.4])
lignin_froot = np.array([0.25,	0.2,	0.2,	0.2,	0.2,	0.1,	0.1,	0.1,	0.1,	0.1,	0.25,	0.2,	0.2,	0.15,	0.15,	0.25,	0.1])


iveg=dataset_static['iveg'][:]

mask = (~np.isnan(iveg)) & (iveg >= 1) & (iveg <= 17)

lignin_leaf_grid  = np.full_like(iveg, np.nan, dtype=float)
lignin_CWD_grid   = np.full_like(iveg, np.nan, dtype=float)
lignin_froot_grid = np.full_like(iveg, np.nan, dtype=float)

lignin_leaf_grid[mask]  =lignin_leaf[iveg[mask].astype(int)-1]
lignin_CWD_grid[mask]   =lignin_CWD[iveg[mask].astype(int)-1]
lignin_froot_grid[mask] =lignin_froot[iveg[mask].astype(int)-1]


new = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/lignin_C_cable.nc", 'w')

new.createDimension("lat", 112 )
new.createDimension("lon", 192 )
new.createDimension("patch", 11 )


new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=dataset_static['y'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=dataset_static['x'][:]

new.createVariable("lignin_leaf_grid", "float64", (u'patch',u'lat',u'lon',))
new.variables['lignin_leaf_grid'].__setattr__("long_name","Lignin:C ratio leaf litter fall")
new.variables['lignin_leaf_grid'][:]=lignin_leaf_grid[:]

new.createVariable("lignin_CWD_grid", "float64", (u'patch',u'lat',u'lon',))
new.variables['lignin_CWD_grid'].__setattr__("long_name","Lignin:C ratio non-leaf litter fall")
new.variables['lignin_CWD_grid'][:]=lignin_CWD_grid[:]

new.createVariable("lignin_froot_grid", "float64", (u'patch',u'lat',u'lon',))
new.variables['lignin_froot_grid'].__setattr__("long_name","Lignin:C ratio of belowground litter fall")
new.variables['lignin_froot_grid'][:]=lignin_froot_grid[:]


