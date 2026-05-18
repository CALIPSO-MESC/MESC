import netCDF4 as nc
import pandas as pd
import numpy as np
import xarray as xr

dataset_daily = nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time-varying variables_HWSD_soil_layer.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/cable_out_hwsd90_cn_1990_extracted_var.nc")
dataset_invariant=nc.Dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc")
USDA_soil_data= nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_mask_cable.nc")

#================================load data==================================================#
lat               =np.array(USDA_soil_data['lat'][:])
lon               =np.array(USDA_soil_data['lon'][:])

USDA_SoilSuborder =np.array(USDA_soil_data['USDA_SoilSuborder'][:])
max_PFT           =np.array(USDA_soil_data['max_PFT'][:])
max_PFTfrac       =np.array(USDA_soil_data['max_PFTfrac'][:])

# npp(PFT, lat, lon)
npp         =np.array(dataset_invariant['npp'][:])
HWSD_SOC    =np.array(dataset_invariant['HWSD_SOC'][:])
HWSD_bulk_density=np.array(dataset_invariant['HWSD_bulk_density'][:])
PFTfrac     =np.array(dataset_invariant['PFTfrac'][:])

Ald         =np.array(dataset_invariant['Ald'][:])
Alo         =np.array(dataset_invariant['Alo'][:])
Fed         =np.array(dataset_invariant['Fed'][:])
Feo         =np.array(dataset_invariant['Feo'][:])
Ca          =np.array(dataset_invariant['Ca'][:])
Mg          =np.array(dataset_invariant['Mg'][:])

pH          =np.array(dataset_invariant['pH'][:])
area        =np.array(dataset_invariant['area'][:])
rhosoil     =np.array(dataset_invariant['rhosoil'][:])

bch         =np.array(dataset_static['bch'][0,:])
clay        =np.array(dataset_static['clay'][0,:])
isoil       =np.array(dataset_static['isoil'][0,:])
sand        =np.array(dataset_static['sand'][0,:])
sfc         =np.array(dataset_static['sfc'][0,:])
silt        =np.array(dataset_static['silt'][0,:])
ssat        =np.array(dataset_static['ssat'][0,:])
sucs        =np.array(dataset_static['sucs'][0,:])
swilt       =np.array(dataset_static['swilt'][0,:])



# SoilTemp(time, soil_hwsd, PFT, lat, lon)
SoilTemp  =np.array(dataset_daily['SoilTemp'][:])
SoilMoist =np.array(dataset_daily['SoilMoist'][:])

cleaf2met =np.array(dataset_daily['cleaf2met'][:])
cleaf2str =np.array(dataset_daily['cleaf2str'][:])
croot2met =np.array(dataset_daily['croot2met'][:])
croot2str =np.array(dataset_daily['croot2str'][:])
cwood2cwd =np.array(dataset_daily['cwood2cwd'][:])

nleaf2met =np.array(dataset_daily['nleaf2met'][:])
nleaf2str =np.array(dataset_daily['nleaf2str'][:])
nroot2met =np.array(dataset_daily['nroot2met'][:])
nroot2str =np.array(dataset_daily['nroot2str'][:])
nwood2cwd =np.array(dataset_daily['nwood2cwd'][:])

Belowground_litter_fall        =np.array(dataset_daily['Belowground_litter_fall'][:])
Leaf_fall                      =np.array(dataset_daily['Leaf_fall'][:])
non_leaf_aboveground_litterfall=np.array(dataset_daily['non_leaf_aboveground_litterfall'][:])

water_potential=np.array(dataset_daily['water_potential'][:])

PFTfrac[PFTfrac == -99999]     = np.nan
SoilTemp[SoilTemp <= -99999]   = np.nan
SoilMoist[SoilMoist == -99999] = np.nan

cleaf2met[cleaf2met == -99999] = np.nan
cleaf2str[cleaf2str == -99999] = np.nan
croot2met[croot2met == -99999] = np.nan
croot2str[croot2str == -99999] = np.nan
cwood2cwd[cwood2cwd == -99999] = np.nan

nleaf2met[nleaf2met == -99999] = np.nan
nleaf2str[nleaf2str == -99999] = np.nan
nroot2met[nroot2met == -99999] = np.nan
nroot2str[nroot2str == -99999] = np.nan
nwood2cwd[nwood2cwd == -99999] = np.nan

Belowground_litter_fall[Belowground_litter_fall == -99999] = np.nan
Leaf_fall[Leaf_fall == -99999] = np.nan
non_leaf_aboveground_litterfall[non_leaf_aboveground_litterfall == -99999] = np.nan

npp[npp == -99999] = np.nan
HWSD_SOC[HWSD_SOC == -99999] = np.nan
water_potential[water_potential <= -9999] = np.nan


npp=np.nansum(npp*PFTfrac[:],axis=0)
SoilTemp=np.nanmean(SoilTemp*PFTfrac[np.newaxis,np.newaxis,:],axis=2)
SoilMoist=np.nanmean(SoilMoist*PFTfrac[np.newaxis,np.newaxis,:],axis=2)

cleaf2met=np.nansum(cleaf2met*PFTfrac[np.newaxis,:],axis=1)
cleaf2str=np.nansum(cleaf2str*PFTfrac[np.newaxis,:],axis=1)
croot2met=np.nansum(croot2met*PFTfrac[np.newaxis,:],axis=1)
croot2str=np.nansum(croot2str*PFTfrac[np.newaxis,:],axis=1)
cwood2cwd=np.nansum(cwood2cwd*PFTfrac[np.newaxis,:],axis=1)

nleaf2met=np.nansum(nleaf2met*PFTfrac[np.newaxis,:],axis=1)
nleaf2str=np.nansum(nleaf2str*PFTfrac[np.newaxis,:],axis=1)
nroot2met=np.nansum(nroot2met*PFTfrac[np.newaxis,:],axis=1)
nroot2str=np.nansum(nroot2str*PFTfrac[np.newaxis,:],axis=1)
nwood2cwd=np.nansum(nwood2cwd*PFTfrac[np.newaxis,:],axis=1)

Belowground_litter_fall=np.nansum(Belowground_litter_fall*PFTfrac[np.newaxis,:,:,:],axis=1)
Leaf_fall=np.nansum(Leaf_fall*PFTfrac[np.newaxis,:,:,:],axis=1)
non_leaf_aboveground_litterfall=np.nansum(non_leaf_aboveground_litterfall*PFTfrac[np.newaxis,:,:,:],axis=1)

#================================extra data==================================================#
soil_mask = (USDA_SoilSuborder != -1) & (pH>3) & (npp<1500)

filtered_soil_types = USDA_SoilSuborder[soil_mask]
unique_soil_types, counts = np.unique(filtered_soil_types, return_counts=True)

for soil_type, count in zip(unique_soil_types, counts):
    print(f"土壤类型 {soil_type}: {count} 个像元")

print("满足条件的像元总数:"+str(np.sum(soil_mask)) )
    
filtered_PFT_types = max_PFT[soil_mask]
unique_PFT_types, counts = np.unique(filtered_PFT_types, return_counts=True)

for PFT_type, count in zip(unique_PFT_types, counts):
    print(f"PFT类型 {PFT_type}: {count} 个像元")



lat_indices, lon_indices = np.where(soil_mask)

lat_flat=lat[lat_indices]
lon_flat=lon[lon_indices]

max_PFT_flat=max_PFT[soil_mask]
max_PFTfrac_flat=max_PFTfrac[soil_mask]
USDA_SoilSuborder_flat=USDA_SoilSuborder[soil_mask]
print(npp.shape)
npp_flat=npp[soil_mask]
HWSD_SOC_flat=HWSD_SOC[:,soil_mask]
HWSD_bulk_density_flat=HWSD_bulk_density[:,soil_mask]

Ald_flat=Ald[soil_mask]
Alo_flat=Alo[soil_mask]
Fed_flat=Fed[soil_mask]
Feo_flat=Feo[soil_mask]
Ca_flat=Ca[soil_mask]
Mg_flat=Mg[soil_mask]
pH_flat =pH[soil_mask]
area_flat=area[soil_mask]

SoilTemp_flat =SoilTemp[:, :, soil_mask]
SoilMoist_flat=SoilMoist[:, :, soil_mask]

cleaf2met_flat=cleaf2met[:,soil_mask]
cleaf2str_flat=cleaf2str[:,soil_mask]
croot2met_flat=croot2met[:,soil_mask]
croot2str_flat=croot2str[:,soil_mask]
cwood2cwd_flat=cwood2cwd[:,soil_mask]

nleaf2met_flat=nleaf2met[:,soil_mask]
nleaf2str_flat=nleaf2str[:,soil_mask]
nroot2met_flat=nroot2met[:,soil_mask]
nroot2str_flat=nroot2str[:,soil_mask]
nwood2cwd_flat=nwood2cwd[:,soil_mask]

Belowground_litter_fall_flat=Belowground_litter_fall[:,soil_mask]
Leaf_fall_flat=Leaf_fall[:,soil_mask]
non_leaf_aboveground_litterfall_flat=non_leaf_aboveground_litterfall[:,soil_mask]

water_potential_flat=water_potential[:,:,soil_mask]

#print(soil_mask.shape)
#print(SoilTemp_flat.shape)
#print(SoilMoist_flat.shape)
#print(lat_flat.shape)
#print(lon_flat.shape)
#print(npp.shape)
#print(npp_flat.shape)

bch_flat= bch[soil_mask]
clay_flat= clay[soil_mask]
isoil_flat= isoil[soil_mask]
sand_flat= sand[soil_mask]
sfc_flat= sfc[soil_mask]
silt_flat= silt[soil_mask]
ssat_flat= ssat[soil_mask]
sucs_flat= sucs[soil_mask]
swilt_flat= swilt[soil_mask]
rhosoil_flat= rhosoil[soil_mask]
print(swilt_flat.shape)


#===================================write nc file========================================#

time = dataset_daily['time'][:]/86400+0.5   #seconds to day

soil_hwsd=np.arange(1, 8)

print("write nc file")

ds = xr.Dataset(
    {
        "lat": (["nsite"], lat_flat),
        "lon": (["nsite"], lon_flat),
        "area": (["nsite"], area_flat),
        
        "max_PFT": (["nsite"], max_PFT_flat+1),
        "max_PFTfrac": (["nsite"], max_PFTfrac_flat),
        "USDA_SoilSuborder": (["nsite"], USDA_SoilSuborder_flat),
        "npp": (["nsite"], npp_flat),
        
        "Ald": (["nsite"], Ald_flat),
        "Alo": (["nsite"], Alo_flat),
        "Fed": (["nsite"], Fed_flat),
        "Feo": (["nsite"], Feo_flat),
        "Ca": (["nsite"],  Ca_flat),
        "Mg": (["nsite"],  Mg_flat),

        "pH": (["nsite"],  pH_flat),
       
        "bch": (["nsite"], bch_flat),
        "clay": (["nsite"], clay_flat),
        "isoil": (["nsite"], isoil_flat),
        "sand": (["nsite"], sand_flat),
        "sfc": (["nsite"], sfc_flat),
        "silt": (["nsite"], silt_flat),
        "ssat": (["nsite"], ssat_flat),
        "sucs": (["nsite"], sucs_flat),
        "swilt": (["nsite"], swilt_flat),
        "bulk_density": (["nsite"], rhosoil_flat),

        "HWSD_SOC": (["soil_hwsd", "nsite"], HWSD_SOC_flat),
        "HWSD_bulk_density":(["soil_hwsd", "nsite"], HWSD_bulk_density_flat),
        "SoilMoist": (["time", "soil_hwsd", "nsite"], SoilMoist_flat),
        "SoilTemp": (["time", "soil_hwsd", "nsite"], SoilTemp_flat),
        "water_potential": (["time", "soil_hwsd", "nsite"], water_potential_flat),

        "Belowground_litter_fall": (["time", "nsite"], Belowground_litter_fall_flat),
        "Leaf_fall": (["time", "nsite"], Leaf_fall_flat),
        "non_leaf_aboveground_litterfall": (["time", "nsite"], non_leaf_aboveground_litterfall_flat),
        
        "cleaf2met": (["time", "nsite"], cleaf2met_flat), 
        "cleaf2str": (["time", "nsite"], cleaf2str_flat),
        "croot2met": (["time", "nsite"], croot2met_flat),
        "croot2str": (["time", "nsite"], croot2str_flat),
        "cwood2cwd": (["time", "nsite"], cwood2cwd_flat), 
        
        "nleaf2met": (["time", "nsite"], nleaf2met_flat), 
        "nleaf2str": (["time", "nsite"], nleaf2str_flat),
        "nroot2met": (["time", "nsite"], nroot2met_flat),
        "nroot2str": (["time", "nsite"], nroot2str_flat),
        "nwood2cwd": (["time", "nsite"], nwood2cwd_flat),               

    },
    coords={
        "time": time,
        "soil_hwsd": soil_hwsd,
    },
)

for var in ds:
    ds[var].encoding["_FillValue"] = -99999.00

ds['time'].attrs['long_name'] = "day of year"
ds['soil_hwsd'].attrs['long_name'] = "HWSD soil layer"
ds['soil_hwsd'].attrs['layer1'] = "0-20cm"
ds['soil_hwsd'].attrs['layer2'] = "20-40cm"
ds['soil_hwsd'].attrs['layer3'] = "40-60cm"
ds['soil_hwsd'].attrs['layer4'] = "60-80cm"
ds['soil_hwsd'].attrs['layer5'] = "80-100cm"
ds['soil_hwsd'].attrs['layer6'] = "100-150cm"
ds['soil_hwsd'].attrs['layer7'] = "150-200cm"

ds['lat'].attrs['long_name'] = "latitude"
ds['lon'].attrs['long_name'] = "longitude"
ds['max_PFT'].attrs['long_name'] = "cable max_PFTfrac > 0.5"
ds['max_PFTfrac'].attrs['long_name'] = "cable PFT, when max_PFTfrac > 0.5"

ds['area'].attrs['long_name'] = "land area"
ds['area'].attrs['units'] = "m^2"

ds['bulk_density'].attrs['long_name'] = "bulk density"
ds['bulk_density'].attrs['units'] = "kg/m^3"

ds['Ald'].attrs['long_name'] = "dithionite-extractable aluminum"
ds['Ald'].attrs['units'] = "kg m^-2"
ds['Ald'].attrs['depth'] = "0-20cm"
ds['Alo'].attrs['long_name'] = "oxalate-extractable aluminum"
ds['Alo'].attrs['units'] = "kg m^-2"
ds['Alo'].attrs['depth'] = "0-20cm"
ds['Fed'].attrs['long_name'] = "dithionite-extractable iron"
ds['Fed'].attrs['units'] = "kg m^-2"
ds['Fed'].attrs['depth'] = "0-20cm"
ds['Feo'].attrs['long_name'] = "oxalate-extractable iron"
ds['Feo'].attrs['units'] = "kg m^-2"
ds['Feo'].attrs['depth'] = "0-20cm"
ds['Ca'] .attrs['long_name'] = "Exchangeable calcium"
ds['Ca'] .attrs['units'] = "cmol/kg"
ds['Ca'] .attrs['depth'] = "0-20cm"
ds['Mg'] .attrs['long_name'] = "Exchangeable magnesium"
ds['Mg'] .attrs['units'] = "cmol/kg"
ds['Mg'] .attrs['depth'] = "0-20cm"


ds['pH'].attrs['long_name'] = "soil pH"

ds['SoilTemp'].attrs['long_name'] = "Average layer soil temperature"
ds['SoilTemp'].attrs['units'] = "Celsius degree"

ds['SoilMoist'].attrs['long_name'] = "Average layer soil moisture"
ds['SoilMoist'].attrs['units'] = "m^3/m^3"

ds['water_potential'].attrs['long_name'] = "Soil water potential by HWSD layer"
ds['water_potential'].attrs['units'] = "kPa"

ds['bch'].attrs['long_name'] = "Parameter b, Campbell eqn 1985"
ds['clay'].attrs['long_name'] = "Fraction of soil which is clay"
ds['isoil'].attrs['long_name'] = "Soil type"
ds['sand'].attrs['long_name'] = "Fraction of soil which is sand"
ds['sfc'].attrs['long_name'] = "Fraction of soil volume which is water @ field capacity"
ds['silt'].attrs['long_name'] = "Fraction of soil which is silt"
ds['ssat'].attrs['long_name'] = "Fraction of soil volume which is water @ saturation"
ds['sucs'].attrs['long_name'] = "Suction @ saturation"
ds['swilt'].attrs['long_name'] = "Fraction of soil volume which is water @ wilting point"


ds['npp'].attrs['long_name'] = "mean annual NPP"
ds['npp'].attrs['units'] = "g C/m2/year"
ds['HWSD_SOC'].attrs['long_name'] = "HWSD SOC with HWSD soil layer in CABLE resolution"
ds['HWSD_SOC'].attrs['units'] = "g C/kg soil"
ds['HWSD_bulk_density'].attrs['long_name'] = "HWSD bulk density"
ds['HWSD_bulk_density'].attrs['units'] = "kg/m3"


ds['Belowground_litter_fall'].attrs['long_name'] = "croot2met + croot2str"
ds['Leaf_fall'].attrs['long_name'] = "cleaf2met + cleaf2str"
ds['non_leaf_aboveground_litterfall'].attrs['long_name'] = "cwood2cwd"

#sum npp
sum_leaf = ds['Leaf_fall'].sum(dim='time')  # 形状变为 (nsite,)
sum_non_leaf = ds['non_leaf_aboveground_litterfall'].sum(dim='time')
sum_below = ds['Belowground_litter_fall'].sum(dim='time')
ds['sum_leaf_wood_root_npp'] = sum_leaf + sum_non_leaf + sum_below
ds= ds.where(ds['sum_leaf_wood_root_npp'] > 100, drop=True)
ds['sum_leaf_wood_root_npp'].attrs['long_name'] = "sum leaf_fall, aboveground_litterfall and belowground_litter_fall"
ds['sum_leaf_wood_root_npp'].attrs['units'] = "g C/m2/year"

print(ds)


ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/data_extraction_cable_model_calibration.nc")












