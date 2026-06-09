# Pre-processing for MESC
Provider: Xinyu Zou (zouxinyu@igsnrr.ac.cn)
Date: 2026.05

These scripts are primarily used for processing input data for the MESC model. The MESC model input data include output data from cable and ORCHIDEE, SOC data from HWSD v2, pH data from soilgrids, and iron and aluminum oxide data from ren-gcb2024. These data are eventually reorganized into time-varying and time-invariant datasets. The main processing steps include data format conversion (tif to netcdf), spatial resolution resampling, multi-year averaging, soil depth-weighted conversion based on HWSD soil depth, variable unit conversion, water potential calculation, data integration (merging multiple variables from different sources into a single netcdf file), and value filtering (e.g., PFT frac > 0.5).

Most of the processing is performed using Python scripts, while data format conversion and spatial resolution resampling are performed using GDAL commands. Some variable merging tasks are completed using CDO and NCO commands.
The two folders, cable and ORCHIDEE, correspond to the data processing for each respective model.


## Code Structure 
```
pre-processing/
    |——convert_scale.bash                                       # Resampling and format conversion (tif to netcdf) of HWSD soil property data
    |——available_USDASoilClass_mask.py                          # Generate filter masks for cable and ORCHIDEE
    |——resample_USDA_SoilSuborder.py                            # Resample USDA SoilSuborder    
    |——cable/
        |——add_CABLE_HWSD_SOC.py                                # Add HWSD SOC variable to CABLE_Time_invariant_variables file；
        |——add_CABLE_HWSD_SOC_to_cluster.py                     # Add HWSD SOC variable to calibration_cluster file；
        |——average_PFT.py                                       # Merge all variables by weighted average based on PFT fraction；
        |——cal_lignin_c.py                                      # Map lignin values from lookup table to spatial grid by PFT；
        |——convert_cable_soil_layer_to_HWSD.py                  # Convert CABLE soil layers to HWSD soil layer depths；
        |——convert_npp_unit.py                                  # Convert NPP units；
        |——convert_patch_to_PFT.py                              # Reshape CABLE patch dimension to PFT dimension；
        |——convert_patch_to_PFT_lignin.py                       # Reshape CABLE patch-dimension lignin to PFT dimension；
        |——convert_pnt_to_PFT_CN_ratio.py                       # Reshape CABLE pnt-dimension CN_ratio to PFT dimension；
        |——convert_pnt_to_PFT_litter_fall.py                    # Reshape CABLE pnt-dimension litter-related variables to PFT dimension；
        |——convert_pnt_to_PFT_litter_leaf_root_wood.py          # Reshape CABLE pnt-dimension litter-related variables to PFT dimension(other)；
        |——convert_pnt_to_PFT_npp_Fluxtolitter.py               # Reshape CABLE pnt-dimension NPP-related variables to PFT dimension；
        |——data_extraction_cable_model_calibration.py           # Extract grid-cell variables to point samples based on filter mask；
        |——data_extraction_cable_model_calibration_check.py     # Check units and value ranges of point-sample variables；
        |——fill_nan.py                                          # Replace NaN values with -99999；
        |——HWSD_soil_cable.py                                   # Replace CABLE original soil properties with HWSD soil properties；
        |——update_Al_Fe_Ca_Mg.py                                # Add iron, aluminum oxides, calcium, and magnesium；
        |——water_potential_cal_USDA_Texture_class.py            # Calculate water potential；
    
    |——ORCHIDEE/
        |——add_bulk_density.py                                  # Add bulk density derived from USDA Soil Texture Class lookup table；
        |——add_HWSD_SOC_to_cluster.py                           # Add HWSD SOC variable to calibration_cluster file；
        |——add_litter.py                                        # Add litter variables；
        |——add_npp_var.py                                       # Add NPP variables；
        |——add_ORCHIDEE_HWSD_SOC.py                             # Add HWSD SOC variable to ORCHIDEE_time_invariant file；
        |——add_soilgrid_soil_properties.py                      # Add SoilGrids soil properties variables；
        |——average_PFT.py                                       # Merge all variables by weighted average based on PFT fraction；
        |——convert_npp_unit.py                                  # Convert NPP units；
        |——convert_ORCHIDEE_soil_layer_to_HWSD.py               # Convert ORCHIDEE soil layers to HWSD soil layer depths；
        |——data_extraction_ORCHIDEE_model_calibration.py        # Extract grid-cell variables to point samples based on filter mask；
        |——data_extraction_ORCHIDEE_model_calibration_check.py  # Check units and value ranges of point-sample variables；
        |——fill_nan.py                                          # Replace NaN values with -99999；
        |——HWSD_soil_ORCHIDEE.py                                # Replace ORHCIDEE original soil properties with HWSD soil properties；
        |——update_Al_Fe_Ca_Mg.py                                # Add iron, aluminum oxides, calcium, and magnesium；
        |——water_potential_cal_USDA_Texture_class_ORCHIDEE.py   # Calculate water potential；







