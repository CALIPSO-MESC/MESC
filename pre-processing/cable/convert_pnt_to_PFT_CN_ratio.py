import numpy as np
import netCDF4 as nc


dataset=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/casa_dump_hwsd90_cn_2991-3011_ave_daily.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/casa_out_hwsd90_cn_2990_daily.nc")
dataset_time=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/hwsd90_cn_2990-3009_ave_daily_SoilTemp_SoilMoist.nc")

# CN_ratio_leaf(time, pnt)
CN_ratio_leaf=dataset['CN_ratio_leaf'][:]
CN_ratio_noleaf=dataset['CN_ratio_noleaf'][:]
CN_ratio_belowground=dataset['CN_ratio_belowground'][:]

patchfrac=dataset_static['patchfrac'][:]
iveg=dataset_static['iveg'][:]
time=dataset_time['time'][:]
print(time.shape)





CN_ratio_leaf_new  = np.full((365,17,112,192), np.nan, dtype=float)
CN_ratio_noleaf_new  = np.full((365,17,112,192), np.nan, dtype=float)
CN_ratio_belowground_new= np.full((365,17,112,192), np.nan, dtype=float)

index_lat=(dataset_static['latitude'][:]+55.0)/1.25
index_lon=(dataset_static['longitude'][:]+179.0625)/1.875


for day in range(0,365):
    print(day)
    CN_ratio_leaf_tmp=dataset['CN_ratio_leaf'][day,:]
    CN_ratio_noleaf_tmp=dataset['CN_ratio_noleaf'][day,:]
    CN_ratio_belowground_tmp=dataset['CN_ratio_belowground'][day,:]
    
    data = np.column_stack((iveg,index_lat,index_lon,CN_ratio_leaf_tmp,CN_ratio_noleaf_tmp,CN_ratio_belowground_tmp))
    for iveg_index, y, x, leaf_value,noleaf_value,belowground_value in data:
        CN_ratio_leaf_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = leaf_value
        CN_ratio_noleaf_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = noleaf_value
        CN_ratio_belowground_new[day,iveg_index.astype(int)-1, y.astype(int), x.astype(int)] = belowground_value
    
CN_ratio_leaf_new=np.nanmean(CN_ratio_leaf_new,axis=0)
CN_ratio_noleaf_new=np.nanmean(CN_ratio_noleaf_new,axis=0)
CN_ratio_belowground_new=np.nanmean(CN_ratio_belowground_new,axis=0)
    

new = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CN_ratio_PFT.nc", 'a')


new.createVariable("CN_ratio_leaf", "float64", (u'PFT',u'lat',u'lon',))
new.variables['CN_ratio_leaf'].__setattr__("long_name","C:N ratio leaf litter fall")
new.variables['CN_ratio_leaf'][:]=CN_ratio_leaf_new[:]

new.createVariable("CN_ratio_noleaf", "float64", (u'PFT',u'lat',u'lon',))
new.variables['CN_ratio_noleaf'].__setattr__("long_name","C:N ratio non-leaf litter fall")
new.variables['CN_ratio_noleaf'][:]=CN_ratio_noleaf_new[:]

new.createVariable("CN_ratio_belowground", "float64", (u'PFT',u'lat',u'lon',))
new.variables['CN_ratio_belowground'].__setattr__("long_name","C:N ratio of belowground litter fall")
new.variables['CN_ratio_belowground'][:]=CN_ratio_belowground_new[:]

