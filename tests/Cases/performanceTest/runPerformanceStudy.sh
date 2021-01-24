#!/bin/bash

set -e

cd src/ && wmake && cd ..


cd Case/ && blockMesh 

LOGNAME="performanceRun.log"


if [ -d "constant/WENOBase3" ]; then 
 rm -r constant/WENOBase3
fi

if [ -e "${LOGNAME}" ]; then
    rm ${LOGNAME}
fi

if [ -e "plotPerformance.dat" ]; then
    rm plotPerformance.dat
fi

touch ${LOGNAME}
touch plotPerformance.dat

for i in {1..30}; do
    ../src/performanceRun.exe  >> ${LOGNAME}
    # Store data in file to plot
    buildUp=$(grep "Duration Build-Up:" performanceRun.log | tail -n1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
    runTimes=$(grep "Duration Run-Time:" performanceRun.log | tail -n1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?')
    echo -e "${runTimes}\t${buildUp}" >> plotPerformance.dat
done


 
