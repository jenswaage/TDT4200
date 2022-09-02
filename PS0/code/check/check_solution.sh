echo
echo "Comparing serial solution with reference data"
echo "------------------------------------------------------------"
echo "Grid size: 1024, Max iteration: 100000, Snapshot frequency: 1000"
echo
./serial 1>/dev/null
./check/compare_solutions 1024 data/00050.bin check/references/n1024/00050.bin
echo
