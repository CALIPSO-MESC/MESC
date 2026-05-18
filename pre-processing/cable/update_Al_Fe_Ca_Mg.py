import numpy as np
import netCDF4 as nc
import xarray as xr

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

dataset_Ald=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/Ald_0_20_cable.nc", decode_times=False)
dataset_Alo=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/Alo_0_20_cable.nc", decode_times=False)
dataset_Fed=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/Fed_0_20_cable.nc", decode_times=False)
dataset_Feo=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/Feo_0_20_cable.nc", decode_times=False)

Ald=dataset_Ald['Band1'].values
Alo=dataset_Alo['Band1'].values
Fed=dataset_Fed['Band1'].values
Feo=dataset_Feo['Band1'].values

#convert unit kg/m2 to cmol/kg
h=0.2
MFe=0.05585  #molar mass of Fe, kg/mol
MAl=0.02698  #molar mass of Al, kg/mol
bd = ds['rhosoil'].values
valid_mask = bd!=-99999
Ald[valid_mask] = (Ald[valid_mask] * 100) / (h * bd[valid_mask] * MAl)
Alo[valid_mask] = (Alo[valid_mask] * 100) / (h * bd[valid_mask] * MAl)
Fed[valid_mask] = (Fed[valid_mask] * 100) / (h * bd[valid_mask] * MFe)
Feo[valid_mask] = (Feo[valid_mask] * 100) / (h * bd[valid_mask] * MFe)


Ald[Ald == 0] = -99999
Alo[Alo == 0] = -99999
Fed[Fed == 0] = -99999
Feo[Feo == 0] = -99999


ds['Ald'] = (('lat', 'lon'), Ald)
ds['Alo'] = (('lat', 'lon'), Alo)
ds['Fed'] = (('lat', 'lon'), Fed)
ds['Feo'] = (('lat', 'lon'), Feo)

ds['Ald'].attrs['depth'] = '0-20cm'
ds['Alo'].attrs['depth'] = '0-20cm'
ds['Fed'].attrs['depth'] = '0-20cm'
ds['Feo'].attrs['depth'] = '0-20cm'

ds['Ald'].attrs['units'] = 'cmol/kg'
ds['Alo'].attrs['units'] = 'cmol/kg'
ds['Fed'].attrs['units'] = 'cmol/kg'
ds['Feo'].attrs['units'] = 'cmol/kg'

dataset_Ca=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/EXCA1_cable.nc", decode_times=False)
dataset_Mg=xr.open_dataset("/data1/zxy/SOC_data_calibration/pH_Soil_Al_and_Fe_oxalates/EXMG1_cable.nc", decode_times=False)

Ca=dataset_Ca['EXCA'].fillna(-9999900).values
Mg=dataset_Mg['EXMG'].fillna(-9999900).values
print(type(Ca))

# data var Scale factor=0.01,
Ca_data=(Ca[0,:]*4.5+Ca[0,:]*4.6+Ca[0,:]*7.5+Ca[0,:]*3.4)/20 *0.01 
Mg_data=(Mg[0,:]*4.5+Mg[0,:]*4.6+Mg[0,:]*7.5+Mg[0,:]*3.4)/20 *0.01

Ca_data[Ca_data == np.nan] = -99999
Mg_data[Mg_data == np.nan] = -99999

ds['Ca'] = (('lat', 'lon'), Ca_data)
ds['Mg'] = (('lat', 'lon'), Mg_data)

ds['Ca'].attrs = {'units': 'cmol/kg', 'long_name': 'Exchangeable calcium', 'depth': '0-20cm'}
ds['Mg'].attrs = {'units': 'cmol/kg', 'long_name': 'Exchangeable magnesium', 'depth': '0-20cm'}


ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables_new.nc")

