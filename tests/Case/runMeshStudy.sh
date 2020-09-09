#!/bin/bash

# Create data file
echo "# Mesh Size, linear error, limitedLinear error, WENO error" > PLOT/data.dat


meshSizes=(25 50 100 150)

for cells in ${meshSizes[@]}; do
    sed -i "s/^nCells.*/nCells ${cells};/" system/blockMeshDict
    blockMesh > /dev/null 2>&1 
    ../src/tests.exe [2D] > log

    # Get the mean error and write to data.dat
    linearError=$(grep -m2 "Mean Error Linear" log | tail -n1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
    limitedLinearError=$(grep -m2 "Mean Error LimitedLinear" log | tail -n1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
    WENOError=$(grep -m2 "Mean Error WENO" log | tail -n1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')

    echo -e "${cells}\t${linearError}\t${limitedLinearError}\t${WENOError}" >> PLOT/data.dat
done

rm log

# Plot the results
cd PLOT/ && ./plotResults.sh



