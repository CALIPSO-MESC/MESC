import numpy as np
import pandas as pd
import netCDF4 as nc
from soiltexture import getTexture


SoilMoist_data=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/hwsd90_cn_2990-3009_ave_daily_SoilMoist.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/cable_out_hwsd90_cn_1990_extracted_var.nc")
dataset_grid_access_in=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/grid_access_in.nc")
USDA_Texture_class_data=nc.Dataset("/data1/zxy/SOC_data_calibration/soil water potential/USDA_Texture_class_cable.nc")

USDA_texture_class=USDA_Texture_class_data['USDA_texture_class'][:]

# SoilTemp(time, soil, patch, y, x)
SoilMoist = SoilMoist_data['SoilMoist'][:]
patchfrac=dataset_static['patchfrac'][:]
# isoil(latitude, longitude)
cable_Zobler_soil_texture=dataset_grid_access_in['isoil'][:]

SoilMoist=np.nansum(SoilMoist*patchfrac[np.newaxis,np.newaxis,:],axis=2)


#=================parameter values corresponding to soil texture (USDA system)==================#

USDA_mapping = {
    'clay': 1,
    'clay loam': 2,
    'loam': 3,
    'loamy sand': 4,
    'sand': 5,
    'sandy clay': 6,
    'sandy clay loam': 7,
    'sandy loam': 8,
    'silt': 9,
    'silt clay': 10,
    'silt clay loam': 11,
    'silt loam': 12
}


qr    =       np.array([0.098,0.079,0.061,0.049,0.053,0.117,0.063,0.039,0.05,0.111,0.09,0.065])
qs    =       np.array([0.459,0.442,0.399,0.39,0.375,0.385,0.384,0.387,0.489,0.481,0.482,0.439])
alpha = 10 ** np.array([-1.825,-1.801,-1.954,-1.459,-1.453,-1.476,-1.676,-1.574,-2.182,-1.79,-2.076,-2.296])
n     = 10 ** np.array([0.098,0.151,0.168,0.242,0.502,0.082,0.124,0.161,0.225,0.121,0.182,0.221])
m = 1 - 1/n


#=================parameter values corresponding to soil texture (USDA system)==================#



#mask = (~np.isnan(cable_Zobler_soil_texture)) & (cable_Zobler_soil_texture >= 1) & (cable_Zobler_soil_texture <= 9)
mask = (~np.isnan(USDA_texture_class)) & (USDA_texture_class >= 1) & (USDA_texture_class <= 12)

qr_grid = np.full_like(USDA_texture_class, np.nan, dtype=float)
qs_grid = np.full_like(USDA_texture_class, np.nan, dtype=float)
alpha_grid = np.full_like(USDA_texture_class, np.nan, dtype=float)
n_grid = np.full_like(USDA_texture_class, np.nan, dtype=float)
m_grid = np.full_like(USDA_texture_class, np.nan, dtype=float)


qr_grid[mask]   =qr[USDA_texture_class[mask].astype(int)-1]

qs_grid[mask]   =qs[USDA_texture_class[mask].astype(int)-1]
alpha_grid[mask]=alpha[USDA_texture_class[mask].astype(int)-1]
n_grid[mask]    =n[USDA_texture_class[mask].astype(int)-1]
m_grid[mask]    =m[USDA_texture_class[mask].astype(int)-1]

print(n_grid.shape)

water_potential=np.empty((365, 6, 112, 192))

for i in range(0,365):
    print(i)
    for j in range(0,6):
        water_potential[i,j,:] = -1 / alpha_grid * (((qs_grid - qr_grid) / (SoilMoist[i,j,:] - qr_grid)) ** (1 / m_grid) - 1) ** (1 / n_grid)
        ### unit conversion from cm H2O to kPa    
        water_potential[i,j,:]=water_potential[i,j,:]*0.098
        
        p_mask = (SoilMoist[i,j,:] > qr_grid) & (SoilMoist[i,j,:] < qs_grid)
        water_potential[i,j][~p_mask]=np.nan

water_potential[water_potential < -1500] = -1500


new = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/water_potential_cable_USDA_Texture.nc", 'w')

new.createDimension("time", 365 )
new.createDimension("lat", 112 )
new.createDimension("lon", 192 )
new.createDimension("soil", 6 )

new.createVariable("time","float64", (u'time',))
new.variables['time'].__setattr__("units","seconds since 1990-01-01 00:00:00")
new.variables['time'][:]=SoilMoist_data['time'][:]

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=SoilMoist_data['y'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=SoilMoist_data['x'][:]

new.createVariable("water_potential", "float64", (u'time',u'soil',u'lat',u'lon',))
new.variables['water_potential'].__setattr__("long_name","water_potential(kPa)")
new.variables['water_potential'][:]=water_potential[:]

new.createVariable("SoilMoist", "float64", (u'time',u'soil',u'lat',u'lon',))
new.variables['SoilMoist'].__setattr__("long_name","SoilMoist")
new.variables['SoilMoist'][:]=SoilMoist[:]

new.createVariable("qr_grid", "float64", (u'lat',u'lon',))
new.variables['qr_grid'].__setattr__("long_name","qr_grid")
new.variables['qr_grid'][:]=qr_grid[:]

new.createVariable("qs_grid", "float64", (u'lat',u'lon',))
new.variables['qs_grid'].__setattr__("long_name","qs_grid")
new.variables['qs_grid'][:]=qs_grid[:]

new.createVariable("alpha_grid", "float64", (u'lat',u'lon',))
new.variables['alpha_grid'].__setattr__("long_name","alpha_grid")
new.variables['alpha_grid'][:]=alpha_grid[:]

new.createVariable("n_grid", "float64", (u'lat',u'lon',))
new.variables['n_grid'].__setattr__("long_name","n_grid")
new.variables['n_grid'][:]=n_grid[:]

new.createVariable("m_grid", "float64", (u'lat',u'lon',))
new.variables['m_grid'].__setattr__("long_name","m_grid")
new.variables['m_grid'][:]=m_grid[:]









