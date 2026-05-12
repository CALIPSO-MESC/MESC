import xarray as xr
import numpy as np
import pandas as pd
import netCDF4 as nc


ds_cable = nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_cable.nc")
cable_static = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc")

cable_Ald=cable_static['Ald'][:]
cable_lat_mask=cable_static['lat'][:][:, np.newaxis]< 65
max_PFTfrac_cable = np.nanmax(cable_static['PFTfrac'][:],axis=0)
max_PFT_index_cable = np.nanargmax(cable_static['PFTfrac'][:],axis=0)
npp=np.array(cable_static['npp'][:])
pH=np.array(cable_static['pH'][:])
npp[npp == -99999] = np.nan
npp=np.nansum(npp*max_PFTfrac_cable[:],axis=0)

#mask_cable = max_PFTfrac_cable > 0.5
mask_cable = (max_PFTfrac_cable > 0.5) & (pH>3) & (npp>100) & (npp<1500) & cable_lat_mask & (cable_Ald!=-99999)

max_PFTfrac_cable[~mask_cable]=np.nan
max_PFT_index_cable[~mask_cable]=-1

cable_usda_masked=ds_cable['USDA_soil_map'][:]
cable_usda_masked[~mask_cable]=-1


print(mask_cable.shape)
print(ds_cable['USDA_soil_map'].shape)
print(cable_usda_masked.shape)

cable_unique_values, cable_counts = np.unique(cable_usda_masked, return_counts=True)


cable_result_df = pd.DataFrame({
   'cable_USDA_Value': cable_unique_values,
   'cable_Count': cable_counts
}).sort_values('cable_USDA_Value')
print("cable")
print (cable_result_df)

cable_unique_values, cable_counts = np.unique(ds_cable['USDA_soil_map'][:], return_counts=True)
cable_result_df = pd.DataFrame({
   'cable_USDA_Value': cable_unique_values,
   'cable_Count': cable_counts
}).sort_values('cable_USDA_Value')
print("cable")
print (cable_result_df)


new = nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_mask_cable.nc", 'w')

new.createDimension("lat", len(cable_static['lat']) )
new.createDimension("lon", len(cable_static['lon']) )

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=cable_static['lat'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=cable_static['lon'][:]

new.createVariable("max_PFTfrac", "float64", (u'lat',u'lon',))
new.variables['max_PFTfrac'].__setattr__("long_name","cable max_PFTfrac > 0.5")
new.variables['max_PFTfrac'][:]=max_PFTfrac_cable[:]

new.createVariable("max_PFT", "int32", (u'lat',u'lon'),fill_value=-1)
new.variables['max_PFT'].__setattr__("long_name","cable PFT, when max_PFTfrac > 0.5")
new.variables['max_PFT'][:]=max_PFT_index_cable[:]

new.createVariable("USDA_SoilSuborder", "int32", (u'lat',u'lon'),fill_value=-1)
new.variables['USDA_SoilSuborder'].__setattr__("long_name","SoilSuborder, when max_PFTfrac > 0.5")
new.variables['USDA_SoilSuborder'][:]=cable_usda_masked[:]

new.close()


#==================================ORCHIDEE====================================================#
ds_ORCHIDEE = nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_0.5d.nc")
PFT_ORCHIDEE = nc.Dataset("/data1/zxy/SOC_data_calibration/PFTmap_1990_orchidee.nc")
data_static=nc.Dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc")

ORCHIDEE_lat_mask=data_static['lat'][:][:, np.newaxis]< 65
# maxvegetfrac(time_counter, veget, lat, lon) 
max_PFTfrac_ORCHIDEE = np.nanmax(PFT_ORCHIDEE['maxvegetfrac'][0,:],axis=0)
max_PFT_index_ORCHIDEE = np.nanargmax(PFT_ORCHIDEE['maxvegetfrac'][0,:],axis=0)
pH_ORCHIDEE=data_static['pH'][:]
pH_ORCHIDEE_mask=pH_ORCHIDEE>3
npp_ORCHIDEE_mask=data_static['npp'][:] > 100

# USDA_soil_map(lat, lon)
ORCHIDEE_usda_masked=np.flip(ds_ORCHIDEE['USDA_soil_map'][:], axis=0)
ORCHIDEE_usda_masked = np.ma.masked_where(max_PFTfrac_ORCHIDEE.mask, ORCHIDEE_usda_masked)
Soil_texture=data_static['Soil_texture'][:]
Ald=data_static['Ald'][:]
count = np.sum(Soil_texture != -99999)
print(count)

#mask1=(ORCHIDEE_usda_masked == 1) & (max_PFTfrac_ORCHIDEE > 0.75)  &  (Soil_texture!=-99999) 
#mask2=(ORCHIDEE_usda_masked == 2) & (max_PFTfrac_ORCHIDEE == 1)    &  (Soil_texture!=-99999)
#mask3=(ORCHIDEE_usda_masked == 3) & (max_PFTfrac_ORCHIDEE > 0.6)   &  (Soil_texture!=-99999)
#mask4=(ORCHIDEE_usda_masked == 4) & (max_PFTfrac_ORCHIDEE > 0.75)  &  (Soil_texture!=-99999)
#mask5=(ORCHIDEE_usda_masked == 5) & (max_PFTfrac_ORCHIDEE > 0.75 )   &  (Soil_texture!=-99999)
#mask7=(ORCHIDEE_usda_masked == 7) & (max_PFTfrac_ORCHIDEE > 0.75)  &  (Soil_texture!=-99999)
#mask8=(ORCHIDEE_usda_masked == 8) & (max_PFTfrac_ORCHIDEE > 0.75)  &  (Soil_texture!=-99999)
#mask9=(ORCHIDEE_usda_masked == 9) & (max_PFTfrac_ORCHIDEE > 0.94)  &  (Soil_texture!=-99999)
#mask10=(ORCHIDEE_usda_masked == 10) & (max_PFTfrac_ORCHIDEE > 0.75)&  (Soil_texture!=-99999)
#mask11=(ORCHIDEE_usda_masked == 11) & (max_PFTfrac_ORCHIDEE > 0.75)&  (Soil_texture!=-99999)


mask1=(ORCHIDEE_usda_masked == 1) & (max_PFTfrac_ORCHIDEE > 0.7)    &  (Soil_texture!=-99999)
mask2=(ORCHIDEE_usda_masked == 2) & (max_PFTfrac_ORCHIDEE > 0.5)    &  (Soil_texture!=-99999)
mask3=(ORCHIDEE_usda_masked == 3) & (max_PFTfrac_ORCHIDEE > 0.5)    &  (Soil_texture!=-99999)
mask4=(ORCHIDEE_usda_masked == 4) & (max_PFTfrac_ORCHIDEE > 0.5)    &  (Soil_texture!=-99999)
mask5=(ORCHIDEE_usda_masked == 5) & (max_PFTfrac_ORCHIDEE > 0.5 )   &  (Soil_texture!=-99999)
mask7=(ORCHIDEE_usda_masked == 7) & (max_PFTfrac_ORCHIDEE > 0.7)    &  (Soil_texture!=-99999)
mask8=(ORCHIDEE_usda_masked == 8) & (max_PFTfrac_ORCHIDEE > 0.7)    &  (Soil_texture!=-99999)
mask9=(ORCHIDEE_usda_masked == 9) & (max_PFTfrac_ORCHIDEE > 0.94)   &  (Soil_texture!=-99999)
mask10=(ORCHIDEE_usda_masked == 10) & (max_PFTfrac_ORCHIDEE > 0.8)  &  (Soil_texture!=-99999)
mask11=(ORCHIDEE_usda_masked == 11) & (max_PFTfrac_ORCHIDEE > 0.6)  &  (Soil_texture!=-99999)

mask_ORCHIDEE= (mask1 | mask2 | mask3 | mask4 | mask5 |
                mask7 | mask8 | mask9 | mask10 | mask11 ) & pH_ORCHIDEE_mask & ORCHIDEE_lat_mask & (Ald!=-99999) & npp_ORCHIDEE_mask

max_PFTfrac_ORCHIDEE[~mask_ORCHIDEE]=np.nan
max_PFT_index_ORCHIDEE[~mask_ORCHIDEE]=-1
ORCHIDEE_usda_masked[~mask_ORCHIDEE]=-1


ORCHIDEE_unique_values, ORCHIDEE_counts = np.unique(ORCHIDEE_usda_masked, return_counts=True)

ORCHIDEE_result_df = pd.DataFrame({
   'ORCHIDEE_USDA_Value': ORCHIDEE_unique_values,
   'ORCHIDEE_Count': ORCHIDEE_counts
}).sort_values('ORCHIDEE_USDA_Value')
print("ORCHIDEE")
print (ORCHIDEE_result_df)

ORCHIDEE_unique_values, ORCHIDEE_counts = np.unique(ds_ORCHIDEE['USDA_soil_map'][:], return_counts=True)

ORCHIDEE_result_df = pd.DataFrame({
   'ORCHIDEE_USDA_Value': ORCHIDEE_unique_values,
   'ORCHIDEE_Count': ORCHIDEE_counts
}).sort_values('ORCHIDEE_USDA_Value')
print("ORCHIDEE")
print (ORCHIDEE_result_df)



new = nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_mask_ORCHIDEE.nc", 'w')

new.createDimension("lat", len(PFT_ORCHIDEE['lat']) )
new.createDimension("lon", len(PFT_ORCHIDEE['lon']) )

new.createVariable("lat", "float64", (u'lat',))
new.variables['lat'].__setattr__("long_name","Latitude")
new.variables['lat'][:]=PFT_ORCHIDEE['lat'][:]

new.createVariable("lon", "float64", (u'lon',))
new.variables['lon'].__setattr__("long_name","Longitude")
new.variables['lon'][:]=PFT_ORCHIDEE['lon'][:]

new.createVariable("max_PFTfrac", "float64", (u'lat',u'lon',))
new.variables['max_PFTfrac'].__setattr__("long_name","ORCHIDEE max_PFTfrac")
new.variables['max_PFTfrac'][:]=max_PFTfrac_ORCHIDEE[:]

new.createVariable("max_PFT", "int32", (u'lat',u'lon'),fill_value=-1)
new.variables['max_PFT'].__setattr__("long_name","ORCHIDEE PFT")
new.variables['max_PFT'][:]=max_PFT_index_ORCHIDEE[:]

new.createVariable("USDA_SoilSuborder", "int32", (u'lat',u'lon'),fill_value=-1)
new.variables['USDA_SoilSuborder'].__setattr__("long_name","SoilSuborder")
new.variables['USDA_SoilSuborder'][:]=ORCHIDEE_usda_masked[:]

new.close()






