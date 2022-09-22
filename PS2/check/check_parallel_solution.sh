echo
echo "Comparing parallel solution with reference data"
echo "-----------------------------------------------------------"
echo "Grid size: 256, Max iteration: 5000, Snapshot frequency: 40"
echo
for i in 1 2 4 8
do
    echo "Running with $i processes:"
    mpirun -n $i --oversubscribe ./parallel 1>/dev/null
    ./check/compare_solutions 256 data/00050.bin check/references/n256/00050.bin
    echo
done
