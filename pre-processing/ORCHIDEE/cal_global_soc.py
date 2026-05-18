import numpy as np
import netCDF4 as nc
import xarray as xr


death=[0.2,0.2,0.2,0.2,0.2,0.5,0.5]

#====================================================CABLE=====================================================#
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

bulk_density =np.array(ds['rhosoil'])
cell_area  =np.array(ds['area'])
#HWSD_SOC   =np.array(ds['HWSD_SOC'])

bulk_density[bulk_density<0]=np.nan
cell_area[cell_area<0]=np.nan
#HWSD_SOC[HWSD_SOC==-99999.0]=np.nan

soc_global_all = np.zeros((112, 192))

for i in range(1,8):
    ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_new_cable_D"+str(i)+".nc", decode_times=False)
    HWSD_SOC   =np.array(ds['Band1'])
    ds_bulk_density = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_bulk_cable_D"+str(i)+".nc")
    bulk_density=np.array(ds_bulk_density['Band1'])*1000
    soc_global=HWSD_SOC*bulk_density*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC*bulk_density*death[i-1]/1e3
    print("cable HWSD bulk",np.nansum(soc_global))


soc_global_var = xr.DataArray(
    soc_global_all,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'soc_global',
        'units': 'kg/m2',
    }
) 
ds['soc_global']=soc_global_var 
  
ds.to_netcdf("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_cable_D1-7.nc")


#====================================================CABLE=====================================================#
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

bulk_density =np.array(ds['rhosoil'])
cell_area  =np.array(ds['area'])
#HWSD_SOC   =np.array(ds['HWSD_SOC'])

bulk_density[bulk_density<0]=np.nan
cell_area[cell_area<0]=np.nan
#HWSD_SOC[HWSD_SOC==-99999.0]=np.nan

soc_global_all = np.zeros((112, 192))

for i in range(1,8):
    ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_new_cable_D"+str(i)+".nc", decode_times=False)
    HWSD_SOC   =np.array(ds['Band1'])
    #ds_bulk_density = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_bulk_cable_D"+str(i)+".nc")
    #bulk_density=np.array(ds_bulk_density['Band1'])*1000
    soc_global=HWSD_SOC*bulk_density*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC*bulk_density*death[i-1]/1e3
    print("cable",np.nansum(soc_global))


soc_global_var = xr.DataArray(
    soc_global_all,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'soc_global',
        'units': 'kg/m2',
    }
) 
ds['soc_global']=soc_global_var 
  
ds.to_netcdf("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_cable_D1-7_cable_bulk.nc")



#====================================================CABLE=====================================================#
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables_new.nc", decode_times=False)

bulk_density =np.array(ds['HWSD_bulk_density'])
cell_area  =np.array(ds['area'])
HWSD_SOC   =np.array(ds['HWSD_SOC'])

bulk_density[bulk_density<0]=np.nan
cell_area[cell_area<0]=np.nan
HWSD_SOC[HWSD_SOC==-99999.0]=np.nan

soc_global_all = np.zeros((112, 192))

for i in range(1,8):
    soc_global=HWSD_SOC[i-1,:]*bulk_density[i-1,:]*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC[i-1,:]*bulk_density[i-1,:]*death[i-1]/1e3
    print("CABLE_Time_invariant_variables_new",np.nansum(soc_global))





#====================================================HWSD_SOC_content_0.5d=====================================================#
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)

bulk_density =np.array(ds['bulk_density'])
cell_area  =np.array(ds['cell_area'])
bulk_density=bulk_density[::-1, ...]

bulk_density[bulk_density<0]=0
cell_area[cell_area<0]=0


soc_global_all = np.zeros((360, 720))

for i in range(1,8):
    ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_0.5d_new_D"+str(i)+".nc", decode_times=False)
    ds_bulk_density = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_bulk_0.5d_D"+str(i)+".nc")
    bulk_density=np.array(ds_bulk_density['Band1'])*1000
    HWSD_SOC   =np.array(ds['Band1'])
    soc_global=HWSD_SOC*bulk_density*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC*bulk_density*death[i-1]/1e3
    print("HWSD_SOC_content_0.5d",np.nansum(soc_global))

soc_global_var = xr.DataArray(
    soc_global_all,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'soc_global',
        'units': 'kg/m2',
    }
)

ds['soc_global']=soc_global_var
#
ds.to_netcdf("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_0.5d_D1-7_ave.nc")

 
    
#====================================================HWSD 0.5d ave=====================================================#    
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)

bulk_density =np.array(ds['bulk_density'])
cell_area  =np.array(ds['cell_area'])
bulk_density=bulk_density[::-1, ...]

bulk_density[bulk_density<0]=0
cell_area[cell_area<0]=0

soc_global_all = np.zeros((360, 720))
    
for i in range(1,8):
    ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_0.5d_new_D"+str(i)+".nc", decode_times=False)
    HWSD_SOC   =np.array(ds['Band1'])
    soc_global=HWSD_SOC*bulk_density*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC*bulk_density*death[i-1]/1e3
    print("HWSD ORCHIDEE:",np.nansum(soc_global))  
        
   
soc_global_var = xr.DataArray(
    soc_global_all,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'soc_global',
        'units': 'kg/m2',
    }
) 
ds['soc_global']=soc_global_var 
  
ds.to_netcdf("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_0.5d_D1-7_ave_ORCHIDEE_bulk.nc")


#====================================================ORCHIDEE=====================================================#
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new.nc", decode_times=False)

bulk_density =np.array(ds['HWSD_bulk_density'])
cell_area  =np.array(ds['cell_area'])
HWSD_SOC   =np.array(ds['HWSD_SOC'])

bulk_density[bulk_density<0]=np.nan
cell_area[cell_area<0]=np.nan
HWSD_SOC[HWSD_SOC==-99999.0]=np.nan

soc_global_all = np.zeros((360, 720))

for i in range(1,8):
    soc_global=HWSD_SOC[i-1,:]*bulk_density[i-1,:]*cell_area*death[i-1]/1e15
    soc_global_all+=HWSD_SOC[i-1,:]*bulk_density[i-1,:]*death[i-1]/1e3
    print("ORCHIDEE_time_invariant_new",np.nansum(soc_global))



   
   
   
   
   
   
   
   
   
   
   
   
   
    
    
    