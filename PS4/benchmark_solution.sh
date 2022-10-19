echo
echo "------------------------------------------------------------"
echo "Benchmark small size"
echo "Grid size: 1024, Max iteration: 100000, Snapshot frequency: 1000"
echo "------------------------------------------------------------------"
echo "Serial"
time $(./serial 1>/dev/null)
echo
echo "CUDA"
time $(./parallel 1>/dev/null)
echo
echo "------------------------------------------------------------"
echo "Benchmark medium size"
echo "Grid size: 4096, Max iteration: 100000, Snapshot frequency: 1000"
echo "------------------------------------------------------------------"
echo "Serial"
time $(./serial -n 4096 1>/dev/null)
echo
echo "CUDA"
time $(./parallel -n 4096 1>/dev/null)
echo
echo "------------------------------------------------------------"
echo "Benchmark large size"
echo "Grid size: 16384, Max iteration: 100000, Snapshot frequency: 1000"
echo "------------------------------------------------------------------"
echo "Serial"
time $(./serial -n 16384 1>/dev/null)
echo
echo "CUDA"
time $(./parallel -n 16384 1>/dev/null)
echo
