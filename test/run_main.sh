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
set -e
# --------------------------------------------------
# Set environment variables
# --------------------------------------------------
export OMP_NUM_THREADS=8
# --------------------------------------------------
# Remove old files
# --------------------------------------------------
rm -f -- *.txt mesc.nml
# --------------------------------------------------
# Configure test cases to be run
# --------------------------------------------------
cases=("frc_f3" "hwsd_cable3")
# --------------------------------------------------
# Loop over test cases
# --------------------------------------------------
rm -rf -- output
mkdir -p -- output

for case in "${cases[@]}"; do
  echo "Running test case '${case}' "
  # --------------------------------------------------
  # Copy namelist (parameter) files into run directory
  # --------------------------------------------------
  cp -- "./input/mesc_${case}.nml" mesc.nml
  # --------------------------------------------------
  # Run the test case
  # --------------------------------------------------
  START="$(date +%s)"
  ./main >"output/outval_${case}.txt"
  DURATION="$(($(date +%s) - START))"
  echo "Time taken: ${DURATION} seconds"

  mv -- fort.91 "output/valsoc_91_${case}.txt"
  mv -- fort.92 "output/valsoc_92_${case}.txt"

  # Run diff, capture return code without triggering set -e exit
  diff -- "benchmark/valsoc_91_${case}.txt" "output/valsoc_91_${case}.txt" >"output/diff_valsoc_91_${case}.txt" || true
  diff -- "benchmark/valsoc_92_${case}.txt" "output/valsoc_92_${case}.txt" >"output/diff_valsoc_92_${case}.txt" || true
done

# --------------------------------------------------
# Report test statuses
# --------------------------------------------------
for case in "${cases[@]}"; do
  pass=true
  for id in 91 92; do
    output_file="output/valsoc_${id}_${case}.txt"
    diff_file="output/diff_valsoc_${id}_${case}.txt"
    # Check the output file exists
    if [ ! -e "${output_file}" ]; then
      pass=false
      break
    fi
    # Check file is not empty
    if [ ! -s "${output_file}" ]; then
      pass=false
      break
    fi
    # Check the diff file exists and is empty
    if [ -s "${diff_file}" ]; then
      pass=false
      break
    fi
  done

  if [ "${pass}" = false ]; then
    echo "FAIL: test case '${case}'!"
  else
    echo "PASS: test case '${case}'!"
  fi
done

# rm output/diff_*.txt
echo "===== Job finished of frc and hwsd: $(date) ====="

# --------------------------------------------------
# Run the global_orchidee case
# --------------------------------------------------
rm -f -- *.txt mesc.nml
case_name="global_orchidee"
echo "Running test case '${case_name}' "

cp -- "./input/mesc_${case_name}.nml" mesc.nml
START="$(date +%s)"
./main >"outval_${case_name}.txt"
DURATION="$(($(date +%s) - START))"
echo "Time taken: ${DURATION} seconds"

cp -- vmic_output.nc  "output/vmic_output_${case_name}.nc"
cp -- miccpool_out.nc "output/miccpool_out_${case_name}.nc"
cp -- "outval_${case_name}.txt" output/

echo "===== Job finished of global_orchidee: $(date) ====="