#!/bin/bash

HOME_DIR="."

if [[ $# -eq 1 ]]; then
    HOME_DIR=$1
fi


echo "Current OpenFOAM version is $WM_PROJECT_VERSION."
if [ -e foamVersionThisIsCompiledFor ]; then
    prevVersion=$(<foamVersionThisIsCompiledFor)
    echo "Previously compiled for OpenFOAM ($prevVersion)"
    if [ "$WM_PROJECT_VERSION" != "$prevVersion" ]
    then
        echo "Different than current OpenFOAM ($WM_PROJECT_VERSION) version."
        echo "   Use ./Allwclean before compiling"
        exit 42
    fi
        unset prevVersion
else
        echo "This is a clean install"
        echo $WM_PROJECT_VERSION >foamVersionThisIsCompiledFor
fi

versionFile=foamVersion4weno.H

echo "here $(pwd)"
${HOME_DIR}/versionRules/makeFoamVersionHeader.py $WM_PROJECT_VERSION >$versionFile.tmp


if [ -e $versionFile ]; then                                                    
    nrDiff=$(diff $versionFile.tmp $versionFile | wc -l | tr -d " ")            
    if [[ $nrDiff > 0 ]]; then                                                  
    echo "$versionFile changed"                                                 
    mv $versionFile.tmp $versionFile                                            
    else                                                                        
    # Make sure that not a complete recompilation is triggered                  
    echo "No change to $versionFile"                                            
    rm $versionFile.tmp                                                         
    fi                                                                          
else                                                                            
    echo "No $versionFile. Generating!"                                           
    mv $versionFile.tmp $versionFile                                            
fi                                          
