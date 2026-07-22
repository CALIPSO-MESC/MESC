#!/bin/bash
# Shell script for running the MESC test suite.

# --------------------------------------------------
# Load environment (users can edit if needed)
# --------------------------------------------------
module purge
module load oneapi23u1 netcdf_intel
# module load intel-fc/2020.4.304
# module unload netcdf
# module load netcdf/4.8.1-intel20

# --------------------------------------------------
# Set environment variables
# --------------------------------------------------
export OMP_NUM_THREADS=8

# --------------------------------------------------
# Remove old files
# --------------------------------------------------
rm -f fort.*
rm -f val*.txt
rm -f params1.txt
rm -f params_val.txt
rm -f case.txt

# --------------------------------------------------
# Configure test cases to be run
# --------------------------------------------------
cases=("f3" "cable3")
runs=("frc" "hwsd")

# --------------------------------------------------
# Loop over test cases
# --------------------------------------------------
mkdir -p output
for i in 0 1; do
  case="${cases[${i}]}"
  run="${runs[${i}]}"
  echo "Running test case '${case}', run '${run}'"

  # --------------------------------------------------
  # Copy parameter files
  # --------------------------------------------------
  cp ./input/case_${run}.txt case.txt
  cp ./input/params1_${run}_${case}.txt params1.txt
  cp ./input/params_val_${run}_${case}.txt params_val.txt

  # --------------------------------------------------
  # Run the test case
  # --------------------------------------------------
  ./main >output/outval_${case}_${run}.txt
  mv fort.91 output/valsoc_91_${case}_${run}.txt
  mv fort.92 output/valsoc_92_${case}_${run}.txt
  diff benchmark/valsoc_91_${case}_${run}.txt output/valsoc_91_${case}_${run}.txt >output/diff_valsoc_91_${case}_${run}.txt
  diff benchmark/valsoc_92_${case}_${run}.txt output/valsoc_92_${case}_${run}.txt >output/diff_valsoc_92_${case}_${run}.txt
done

echo "===== Job finished: $(date) ====="
