echo
echo "Comparing parallel solution with reference data"
echo "------------------------------------------------------------"
echo "Checking thread-block local size correctness"
echo "Grid size: 1024, Max iteration: 100000, Snapshot frequency: 1000"
echo "------------------------------------------------------------------"
echo

./parallel -n 1024 1>/dev/null
./check/compare_solutions 1024 data/00050.bin check/references/n1024/00050.bin
echo

echo
echo "------------------------------------------------------------"
echo "Checking global grid size correctness"
echo "Grid size: 2048, Max iteration: 100000, Snapshot frequency: 1000"
echo "------------------------------------------------------------------"
echo

./parallel -n 2048 1>/dev/null
./check/compare_solutions 2048 data/00050.bin check/references/n2048/00050.bin
echo
