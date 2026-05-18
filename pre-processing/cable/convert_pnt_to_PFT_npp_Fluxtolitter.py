import numpy as np
import netCDF4 as nc
import xarray as xr


dataset=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/casa_out_ave_2990-3009_daily_new_var.nc")
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

# FluxCtolitter(time, mlitter, land)
Cnpp=dataset['Cnpp'][0,:]
FluxCtolitter=dataset['FluxCtolitter'][0,:]
FluxNtolitter=dataset['FluxNtolitter'][0,:]

iveg=dataset['iveg'][:]


#Cnpp_new  = np.full((17,112,192), np.nan, dtype=float)
#FluxCtolitter_new  = np.full((3,17,112,192), np.nan, dtype=float)
#FluxNtolitter_new= np.full((3,17,112,192), np.nan, dtype=float)


Cnpp_new = xr.DataArray(
    np.full((17,112,192), -99999.0),
    dims=('PFT', 'lat', 'lon'),
    coords={'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'mean annual NPP',
        'units': 'g C/m2/year',
    }
)

FluxCtolitter_new = xr.DataArray(
    np.full((3,17,112,192), -99999.0),
    dims=('mlitter','PFT', 'lat', 'lon'),
    coords={'mlitter': np.arange(3),'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'amount of plant litter C into litter pool (metabolic, structural and CWD)',
        'units': 'g C/m2/year',
    }
)

FluxNtolitter_new = xr.DataArray(
    np.full((3,17,112,192), -99999.0),
    dims=('mlitter','PFT', 'lat', 'lon'),
    coords={'mlitter': np.arange(3),'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'amount of plant litter N into litter pool (metabolic, structural and CWD)',
        'units': 'g C/m2/year',
    }
)


index_lat=(dataset['latitude'][:]+55.0)/1.25
index_lon=(dataset['longitude'][:]+179.0625)/1.875

print(FluxCtolitter[0,:].shape)
print(iveg.shape)
print(Cnpp.shape)
print(index_lon.shape)
print(Cnpp)
    
data = np.column_stack((iveg,index_lat,index_lon,Cnpp,FluxCtolitter[0,:],FluxCtolitter[1,:],FluxCtolitter[2,:],FluxNtolitter[0,:],FluxNtolitter[1,:],FluxNtolitter[2,:]))
for iveg_index, y, x, npp_value,FluxCtolitter1_value,FluxCtolitter2_value,FluxCtolitter3_value, FluxNtolitter1_value,FluxNtolitter2_value,FluxNtolitter3_value in data:
    Cnpp_new[iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = npp_value
    FluxCtolitter_new[0,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxCtolitter1_value
    FluxCtolitter_new[1,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxCtolitter2_value
    FluxCtolitter_new[2,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxCtolitter3_value
    FluxNtolitter_new[0,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxNtolitter1_value
    FluxNtolitter_new[1,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxNtolitter2_value
    FluxNtolitter_new[2,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = FluxNtolitter3_value
    

ds['Cnpp']=Cnpp_new
ds['FluxCtolitter']=FluxCtolitter_new
ds['FluxNtolitter']=FluxNtolitter_new

ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables_new.nc")

