echo
echo "Comparing parallel solution with reference data"
echo "-----------------------------------------------------------"
echo "Grid size: 256, Max iteration: 5000, Snapshot frequency: 40"
echo
for i in 4
do
    echo "Running with $i processes:"
    ./pthreads 1>/dev/null
    ./check/compare_solutions 256 data/00120.bin check/references/n256/00120.bin
    echo
done
