# --------------------------------------------------
# Load environment (users can edit if needed)
# --------------------------------------------------
module purge
module load oneapi23u1 netcdf_intel
# module load intel-fc/2020.4.304
# module unload netcdf
# module load netcdf/4.8.1-intel20
# --------------------------------------------------
# remove old files 
# --------------------------------------------------
rm -rf fort.*
rm -rf val*.txt 
rm -rf params1.txt
rm -rf params_val.txt
rm -rf case.txt

#### run 1: frc ###
case="f3"
run="frc"
echo ${case} ${run}
# --------------------------------------------------
# copy parameter files 
# --------------------------------------------------
cp ./input/case_${run}.txt case.txt
cp ./input/params1_${run}_${case}.txt params1.txt
cp ./input/params_val_${run}_${case}.txt params_val.txt

if [ ! -d output ]; then
  mkdir output
fi

# --------------------------------------------------
# run 
# --------------------------------------------------
export OMP_NUM_THREADS=8
./main > output/outval_${case}_${run}.txt
mv fort.91 output/valsoc_91_${case}_${run}.txt
mv fort.92 output/valsoc_92_${case}_${run}.txt

diff benchmark/valsoc_91_${case}_${run}.txt output/valsoc_91_${case}_${run}.txt> output/diff_valsoc_91_${case}_${run}.txt
diff benchmark/valsoc_92_${case}_${run}.txt output/valsoc_92_${case}_${run}.txt> output/diff_valsoc_92_${case}_${run}.txt


#### run 2: hwsd ######
case="cable3"
run="hwsd"
echo ${case} ${run}
# --------------------------------------------------
# copy parameter files 
# --------------------------------------------------
cp ./input/case_${run}.txt case.txt
cp ./input/params1_${run}_${case}.txt params1.txt
cp ./input/params_val_${run}_${case}.txt params_val.txt

if [ ! -d output ]; then
  mkdir output
fi

# --------------------------------------------------
# run 
# --------------------------------------------------
export OMP_NUM_THREADS=8
./main > output/outval_${case}_${run}.txt
mv fort.91 output/valsoc_91_${case}_${run}.txt
mv fort.92 output/valsoc_92_${case}_${run}.txt

diff benchmark/valsoc_91_${case}_${run}.txt output/valsoc_91_${case}_${run}.txt> output/diff_valsoc_91_${case}_${run}.txt
diff benchmark/valsoc_92_${case}_${run}.txt output/valsoc_92_${case}_${run}.txt> output/diff_valsoc_92_${case}_${run}.txt

echo "===== Job finished: $(date) ====="
