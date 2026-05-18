import numpy as np
import netCDF4 as nc
import xarray as xr

dataset=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/casa_dump_hwsd90_cn_2991-3011_ave_daily.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/casa_out_hwsd90_cn_2990_daily.nc")
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_HWSD_soil_layer.nc", decode_times=False)


patchfrac=dataset_static['patchfrac'][:]
iveg=dataset_static['iveg'][:]

index_lat=(dataset_static['latitude'][:]+55.0)/1.25
index_lon=(dataset_static['longitude'][:]+179.0625)/1.875


#===================================C leaf root wood=====================================#
cleaf2met_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

cleaf2str_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

croot2met_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

croot2str_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

cwood2cwd_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)


#===================================N leaf root wood=====================================#
nleaf2met_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

nleaf2str_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

nroot2met_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

nroot2str_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

nwood2cwd_new = xr.DataArray(
    np.full((365,17,112,192), -99999.0),
    dims=('time','PFT', 'lat', 'lon'),
    coords={'time': ds['time'],'PFT': ds['PFT'], 'lat': ds['lat'], 'lon': ds['lon']},)

for day in range(0,365):
    print(day)

    cleaf2met_tmp=dataset['cleaf2met'][day,:]
    cleaf2str_tmp=dataset['cleaf2str'][day,:]
    croot2met_tmp=dataset['croot2met'][day,:]
    croot2str_tmp=dataset['croot2str'][day,:]
    cwood2cwd_tmp=dataset['cwood2cwd'][day,:]
    
    nleaf2met_tmp=dataset['nleaf2met'][day,:]
    nleaf2str_tmp=dataset['nleaf2str'][day,:]
    nroot2met_tmp=dataset['nroot2met'][day,:]
    nroot2str_tmp=dataset['nroot2str'][day,:]
    nwood2cwd_tmp=dataset['nwood2cwd'][day,:]
    
    data = np.column_stack((iveg,index_lat,index_lon,
                            cleaf2met_tmp,cleaf2str_tmp,croot2met_tmp,croot2str_tmp,cwood2cwd_tmp,
                            nleaf2met_tmp,nleaf2str_tmp,nroot2met_tmp,nroot2str_tmp,nwood2cwd_tmp))
    print("column_stack")
    for (iveg_index, y, x,
         cleaf2met_value,cleaf2str_value,croot2met_value,croot2str_value,cwood2cwd_value,
         nleaf2met_value,nleaf2str_value,nroot2met_value,nroot2str_value,nwood2cwd_value) in data:
        

        cleaf2met_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = cleaf2met_value
        cleaf2str_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = cleaf2str_value
        croot2met_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = croot2met_value
        croot2str_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = croot2str_value
        cwood2cwd_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = cwood2cwd_value
        
        nleaf2met_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = nleaf2met_value
        nleaf2str_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = nleaf2str_value
        nroot2met_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = nroot2met_value
        nroot2str_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = nroot2str_value
        nwood2cwd_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = nwood2cwd_value
    
    print("integrated three forloop ")
    

ds['cleaf2met']=cleaf2met_new 
ds['cleaf2str']=cleaf2str_new
ds['croot2met']=croot2met_new
ds['croot2str']=croot2str_new
ds['cwood2cwd']=cwood2cwd_new

ds['nleaf2met']=nleaf2met_new 
ds['nleaf2str']=nleaf2str_new
ds['nroot2met']=nroot2met_new
ds['nroot2str']=nroot2str_new
ds['nwood2cwd']=nwood2cwd_new

ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_HWSD_soil_layer_new.nc")

