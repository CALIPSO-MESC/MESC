import netCDF4 as nc
import numpy as np
import os


current_directory = os.path.dirname(os.path.abspath(__file__))


data=nc.Dataset("/data1/share/orchideeshare/ORCHIDEE_INPUT/USDA_SoilSuborder.nc")


USDA_soil_map=data['USDA_soil_map'][:]

lon=data['nav_lon'][0,:]
lat=data['nav_lat'][:,0]

print(lon)
print(lat)
print(len(lon))
print(len(lat))


new = nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder.nc", 'w')

new.createDimension("lat", len(lat) )
new.createDimension("lon", len(lon) )

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=lat

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=lon

new.createVariable("USDA_soil_map", "float64", (u'lat',u'lon',))
new.variables['USDA_soil_map'].__setattr__("long_name","USDA_soil_map")
new.variables['USDA_soil_map'][:]=USDA_soil_map[:]


