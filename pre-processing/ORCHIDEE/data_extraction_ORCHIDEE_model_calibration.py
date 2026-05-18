import netCDF4 as nc
import pandas as pd
import numpy as np
import xarray as xr
import sys
import os

dataset_daily = nc.Dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_ave_1990_2020_time_varing_HWSD_soil_layer_compress.nc")
dataset_static=nc.Dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc")
USDA_soil_data= nc.Dataset("/data1/zxy/SOC_data_calibration/USDA_SoilSuborder_mask_ORCHIDEE.nc")


#================================load data==================================================#
ORCHIDEE_Soil_texture_class=np.array(dataset_static['USDA_Soil_texture_class'])

lat              =np.array(USDA_soil_data['lat'][:])
lon              =np.array(USDA_soil_data['lon'][:])
USDA_SoilSuborder=np.array(USDA_soil_data['USDA_SoilSuborder'][:])
max_PFT          =np.array(USDA_soil_data['max_PFT'][:])
max_PFTfrac      =np.array(USDA_soil_data['max_PFTfrac'][:])

patchfrac        =np.array(dataset_static['maxvegetfrac'][0,:])
npp              =np.array(dataset_static['npp'][:])
HWSD_SOC         =np.array(dataset_static['HWSD_SOC'][:])
HWSD_bulk_density=np.array(dataset_static['HWSD_bulk_density'][:])

Ald              =np.array(dataset_static['Ald'][:])
Alo              =np.array(dataset_static['Alo'][:])
Fed              =np.array(dataset_static['Fed'][:])
Feo              =np.array(dataset_static['Feo'][:])
Ca               =np.array(dataset_static['Ca'][:])
Mg               =np.array(dataset_static['Mg'][:])
pH               =np.array(dataset_static['pH'][:])
cell_area        =np.array(dataset_static['cell_area'][:])

pH_soilgrid               =np.array(dataset_static['pH_soilgrid'][:])
clay_soilgrid               =np.array(dataset_static['clay_soilgrid'][:])
silt_soilgrid               =np.array(dataset_static['silt_soilgrid'][:])


# SoilTemp(time, soil, patch, y, x)
print("load SoilTemp")
sys.stdout.flush()
SoilTemp                       =np.array(dataset_daily['SoilTemp'][:])
print("load SoilMoist")
sys.stdout.flush()
SoilMoist                      =np.array(dataset_daily['SoilMoist'][:])
Belowground_litter_fall        =np.array(dataset_daily['Belowground_litter_fall'][:])
Leaf_fall                      =np.array(dataset_daily['Leaf_fall'][:])
non_leaf_aboveground_litterfall=np.array(dataset_daily['non_leaf_aboveground_litterfall'][:])
water_potential                =np.array(dataset_daily['water_potential'][:])

#print(HWSD_SOC.shape)
#sys.stdout.flush()

soil_mask = USDA_SoilSuborder != -1
lat_indices, lon_indices = np.where(soil_mask)

SoilTemp[SoilTemp == -99999] = np.nan
SoilMoist[SoilMoist == -99999] = np.nan
Belowground_litter_fall[Belowground_litter_fall == -99999] = np.nan
Leaf_fall[Leaf_fall == -99999] = np.nan
non_leaf_aboveground_litterfall[non_leaf_aboveground_litterfall == -99999] = np.nan
water_potential[water_potential <= -1500] = -1500

SoilMoist=np.nansum(SoilMoist*patchfrac[np.newaxis,:,np.newaxis,:,:],axis=1)
Belowground_litter_fall=np.nansum(Belowground_litter_fall*patchfrac[np.newaxis,:,:,:],axis=1)
Leaf_fall=np.nansum(Leaf_fall*patchfrac[np.newaxis,:,:,:],axis=1)
non_leaf_aboveground_litterfall=np.nansum(non_leaf_aboveground_litterfall*patchfrac[np.newaxis,:,:,:],axis=1)



#================================extrat data==================================================#
lat_flat=lat[lat_indices]
lon_flat=lon[lon_indices]
USDA_SoilSuborder_flat=USDA_SoilSuborder[soil_mask]
max_PFT_flat=max_PFT[soil_mask]
max_PFTfrac_flat=max_PFTfrac[soil_mask]

ORCHIDEE_Soil_texture_class_flat=ORCHIDEE_Soil_texture_class[soil_mask]
print("aa")
print(ORCHIDEE_Soil_texture_class_flat)
sys.stdout.flush()

Ald_flat=Ald[soil_mask]
Alo_flat=Alo[soil_mask]
Fed_flat=Fed[soil_mask]
Feo_flat=Feo[soil_mask]
Ca_flat=Ca[soil_mask]
Mg_flat=Mg[soil_mask]
pH_flat =pH[soil_mask]
cell_area_flat=cell_area[soil_mask]

npp_flat=npp[soil_mask]
HWSD_SOC_flat=HWSD_SOC[:,soil_mask]
HWSD_bulk_density_flat=HWSD_bulk_density[:,soil_mask]

SoilTemp_flat =SoilTemp[:, :, soil_mask]
SoilMoist_flat=SoilMoist[:, :, soil_mask]
water_potential_flat=water_potential[:,:,soil_mask]

Belowground_litter_fall_flat=Belowground_litter_fall[:,soil_mask]
Leaf_fall_flat=Leaf_fall[:,soil_mask]
non_leaf_aboveground_litterfall_flat=non_leaf_aboveground_litterfall[:,soil_mask]

pH_soilgrid_flat =pH_soilgrid[soil_mask]
clay_soilgrid_flat =clay_soilgrid[soil_mask]
silt_soilgrid_flat =silt_soilgrid[soil_mask]


#print(HWSD_SOC_flat.shape)
#sys.stdout.flush()
#
#print(soil_mask.shape)
#print(SoilTemp_flat.shape)
#print(SoilMoist_flat.shape)
#print(lat_flat.shape)
#print(lon_flat.shape)

clay=np.array([0.03,	0.06,	0.11,	0.19,	0.1,	0.2,	0.27,	0.33,	0.33,	0.41,	0.46,	0.55])
silt=np.array([0.04,	0.13,	0.26,	0.64,	0.84,	0.4,	0.19,	0.59,	0.37,	0.11,	0.48,	0.3])

soil_porosity=np.array([0.43,	0.41,	0.41,	0.45,	0.46,	0.43,	0.39,	0.43,	0.41,	0.38,	0.36,	0.38])
#bulk density=（1-porosity）* particle density(2.65)
bulk_density=np.array([1.5105,	1.5635,	1.5635,	1.4575,	1.431,	1.5105,	1.6165,	1.5105,	1.5635,	1.643,	1.696,	1.643])*1000
Volumetric_water_content_at_field_capacity=np.array([0.0493,	0.071,	0.1218,	0.2402,	0.2582,	0.1654,	0.1695,	0.3383,	0.2697,	0.2672,	0.337,	0.3469])
Volumetric_water_content_at_wilting_point =np.array([0.045,	0.057,	0.0657,	0.1039,	0.0901,	0.0884,	0.1112,	0.1967,	0.1496,	0.1704,	0.2665,	0.2707])


print(ORCHIDEE_Soil_texture_class_flat)
print(np.max(ORCHIDEE_Soil_texture_class_flat))
print(np.min(ORCHIDEE_Soil_texture_class_flat))

clay_flat=clay[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]
silt_flat=silt[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]
soil_porosity_flat = soil_porosity[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]
bulk_density_flat=bulk_density[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]
sfc_flat=Volumetric_water_content_at_field_capacity[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]
swilt_flat=Volumetric_water_content_at_wilting_point[ORCHIDEE_Soil_texture_class_flat.astype(int) - 1]


#===================================write nc file========================================#

time = dataset_daily['time'][:]/86400+0.5   #seconds to day

soil_hwsd=np.arange(1, 8)

print("write nc file")
sys.stdout.flush()

ds = xr.Dataset(
    {
        "lat": (["nsite"], lat_flat),
        "lon": (["nsite"], lon_flat),
        "cell_area": (["nsite"], cell_area_flat),
        
        "max_PFT": (["nsite"], max_PFT_flat+1),
        "max_PFTfrac": (["nsite"], max_PFTfrac_flat),
        "USDA_SoilSuborder": (["nsite"], USDA_SoilSuborder_flat),
        "USDA_Soil_texture_class": (["nsite"], ORCHIDEE_Soil_texture_class_flat),
        "clay_ORCHIDEE": (["nsite"], clay_flat),
        "silt_ORCHIDEE": (["nsite"], silt_flat),
        "clay_soilgrid": (["nsite"], clay_soilgrid_flat),
        "silt_soilgrid": (["nsite"], silt_soilgrid_flat),

        "Ald": (["nsite"], Ald_flat),
        "Alo": (["nsite"], Alo_flat),
        "Fed": (["nsite"], Fed_flat),
        "Feo": (["nsite"], Feo_flat),
        "Ca": (["nsite"],  Ca_flat),
        "Mg": (["nsite"],  Mg_flat),

        "pH":  (["nsite"], pH_flat),
        "pH_soilgrid":  (["nsite"], pH_soilgrid_flat),
        
        "npp": (["nsite"], npp_flat),

        "soil_porosity_ORCHIDEE": (["nsite"], soil_porosity_flat),
        "bulk_density_ORCHIDEE": (["nsite"], bulk_density_flat),
        "sfc": (["nsite"], sfc_flat),
        "swilt": (["nsite"], swilt_flat),
        

        "HWSD_SOC": (["soil_hwsd", "nsite"], HWSD_SOC_flat),
        "HWSD_bulk_density":(["soil_hwsd", "nsite"], HWSD_bulk_density_flat),
        "SoilMoist": (["time", "soil_hwsd", "nsite"], SoilMoist_flat),
        "SoilTemp": (["time", "soil_hwsd", "nsite"], SoilTemp_flat),
        "water_potential": (["time", "soil_hwsd","nsite"], water_potential_flat),

        "Belowground_litter_fall": (["time", "nsite"], Belowground_litter_fall_flat),
        "Leaf_fall": (["time", "nsite"], Leaf_fall_flat),
        "non_leaf_aboveground_litterfall": (["time", "nsite"], non_leaf_aboveground_litterfall_flat),


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
ds['max_PFT'].attrs['long_name'] = "max Vegetation Types"
ds['max_PFTfrac'].attrs['long_name'] = "Fraction of Vegetation Types"

ds['cell_area'].attrs['long_name'] = "area of grid cell"
ds['cell_area'].attrs['units'] = "m^2"

ds['bulk_density_ORCHIDEE'].attrs['long_name'] = "Soil bulk density"
ds['bulk_density_ORCHIDEE'].attrs['units'] = "kg/m3"
ds['bulk_density_ORCHIDEE'].attrs['description'] = "soil porosity is available，particle density=2.65 g/cm^3, bulk density=（1-porosity）* particle density"

ds['soil_porosity_ORCHIDEE'].attrs['long_name'] = "soil porosity"
ds['soil_porosity_ORCHIDEE'].attrs['description'] = "soil texture-related"

ds['USDA_Soil_texture_class'].attrs['long_name'] = "USDA Soil texture class"
ds['USDA_SoilSuborder'].attrs['long_name'] = "USDA SoilSuborder"
ds['clay_ORCHIDEE'].attrs['long_name'] = "Clay fraction"
ds['clay_ORCHIDEE'].attrs['description'] = "soil texture-related"
ds['silt_ORCHIDEE'].attrs['long_name']   = "Silt fraction"
ds['silt_ORCHIDEE'].attrs['description']   = "soil texture-related"

ds['swilt'].attrs['long_name'] = "soil water content at wilting point"
ds['swilt'].attrs['units'] = "vol water/vol soil"
ds['swilt'].attrs['description'] = "soil texture-related"
ds['sfc'].attrs['long_name'] = "Soil water content at field capacity"
ds['sfc'].attrs['units'] = "vol water/vol soil"
ds['sfc'].attrs['description'] = "soil texture-related"


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

ds['pH_soilgrid'].attrs['depth'] = "0-20cm"
ds['clay_soilgrid'].attrs['depth'] = "0-20cm"
ds['silt_soilgrid'].attrs['depth'] = "0-20cm"

ds['pH_soilgrid'].attrs['units'] = "g/100g(%)"
ds['clay_soilgrid'].attrs['units'] = "g/100g(%)"
ds['silt_soilgrid'].attrs['units'] = "g/100g(%)"

ds['water_potential'].attrs['long_name'] = "Soil water potential by HWSD layer"
ds['water_potential'].attrs['units'] = "kPa"

ds['SoilTemp'].attrs['long_name'] = "Average layer soil temperature"
ds['SoilTemp'].attrs['units'] = "Celsius degree"

ds['SoilMoist'].attrs['long_name'] = "Average layer soil moisture"
ds['SoilMoist'].attrs['units'] = "m^3/m^3"

ds['npp'].attrs['long_name'] = "mean annual NPP"
ds['npp'].attrs['units'] = "g C/m2/year"
ds['HWSD_SOC'].attrs['long_name'] = "HWSD SOC with HWSD soil layer in ORCHIDEE resolution"
ds['HWSD_SOC'].attrs['units'] = "g C/kg soil"
ds['HWSD_bulk_density'].attrs['long_name'] = "HWSD bulk density"
ds['HWSD_bulk_density'].attrs['units'] = "kg/m3"

ds['Belowground_litter_fall'].attrs['long_name'] = "HEART_BE_BM_LITTER + ROOT_BM_LITTER + SAP_BE_BM_LITTER"
ds['Leaf_fall'].attrs['long_name'] = "LEAF_BM_LITTER"
ds['non_leaf_aboveground_litterfall'].attrs['long_name'] = "FRUIT_BM_LITTER + HEART_AB_BM_LITTER + RESERVE_BM_LITTER + SAP_AB_BM_LITTER"
ds['Belowground_litter_fall'].attrs['units'] = "gC/m^2/day"
ds['Leaf_fall'].attrs['units'] = "gC/m^2/day"
ds['non_leaf_aboveground_litterfall'].attrs['units'] = "gC/m^2/day"

#sum npp
sum_leaf = ds['Leaf_fall'].sum(dim='time')  # 形状变为 (nsite,)
sum_non_leaf = ds['non_leaf_aboveground_litterfall'].sum(dim='time')
sum_below = ds['Belowground_litter_fall'].sum(dim='time')
ds['sum_leaf_wood_root_npp'] = sum_leaf + sum_non_leaf + sum_below
ds= ds.where(ds['sum_leaf_wood_root_npp'] > 100, drop=True)
ds['sum_leaf_wood_root_npp'].attrs['long_name'] = "sum leaf_fall, aboveground_litterfall and belowground_litter_fall"
ds['sum_leaf_wood_root_npp'].attrs['units'] = "g C/m2/year"

print(ds)


ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/data_extraction_ORCHIDEE_model_calibration.nc")





