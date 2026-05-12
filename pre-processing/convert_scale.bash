
soilprop="coarse"
for i in $(seq 1 7); do
    gdalwarp -r average -ts 720 360 -of NetCDF "/data1/zxy/SOC_data_calibration/HWSD/HWSD_${soilprop}_D${i}.tif" "/data1/zxy/SOC_data_calibration/HWSD/HWSD_${soilprop}_0.5d_D${i}.nc"
    gdalwarp -r average -ts 192 112 -te -180 -55 180 85 -of NetCDF "/data1/zxy/SOC_data_calibration/HWSD/HWSD_${soilprop}_D${i}.tif" "/data1/zxy/SOC_data_calibration/HWSD/HWSD_${soilprop}_cable_D${i}.nc"
    
done