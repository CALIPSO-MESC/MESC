#!/bin/bash
# Shell script for running the MESC test suite.

# --------------------------------------------------
# Load environment
#
# This will be specific to the machine you are working on and will need to be
# completed by the user. For example:
# ```
# module purge
# module load oneapi23u1   # load intel compiler
# module load netcdf_intel # load netcdf library
# ```
# --------------------------------------------------

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
cases=("frc_f3" "hwsd_cable3" "orchidee_global")

# --------------------------------------------------
# Loop over test cases
# --------------------------------------------------
mkdir -p output
for i in {0..2}; do
  case="${cases[${i}]}"

  echo "Running test case '${case}' "

  # --------------------------------------------------
  # Copy parameter files
  # --------------------------------------------------
  cp ./input/mesc_${case}.nml mesc.nml
  cp ./input/parameters_${case}.txt parameters.txt
  cp ./input/params_val_${case}.txt params_val.txt

  # --------------------------------------------------
  # Run the test case
  # --------------------------------------------------
  ./main >output/outval_${case}.txt
  mv fort.91 output/valsoc_91_${case}.txt
  mv fort.92 output/valsoc_92_${case}.txt
  diff benchmark/valsoc_91_${case}.txt output/valsoc_91_${case}.txt >output/diff_valsoc_91_${case}.txt
  diff benchmark/valsoc_92_${case}.txt output/valsoc_92_${case}.txt >output/diff_valsoc_92_${case}.txt
done

# --------------------------------------------------
# Report test statuses
# --------------------------------------------------
for i in {0..2}; do
  case="${cases[${i}]}"
  pass=1
  for id in 91 92; do
    if [ "$(cat output/diff_valsoc_${id}_${case}.txt)" ]; then
      pass=0
      break
    fi
  done
  if [ ${pass} ]; then
    echo "PASS: test case '${case}'"
  else
    echo "FAIL: test case '${case}'"
  fi
done
rm output/diff_*.txt
echo "===== Job finished: $(date) ====="
