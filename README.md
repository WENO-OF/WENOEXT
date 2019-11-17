WENO framework
====================
Weighted essentially non-oscillatory library for the framework of OpenFOAM.

Tested versions:
    - OpenFOAM 2.3.x
    - OpenFOAM-dev
	

Installation
============
1. rename the folder in "WENOEXT"
2. move the folder "WENOEXT" into $FOAM_SRC
3. Optional: 
      3.1: add the following line to your ~/.bashrc after ". ~/OpenFOAM/OpenFOAM-2.3.x/bashrc":
           . $FOAM_SRC/WENOEXT/bashrc
      3.2: parse your ~/.bashrc or open a new terminal
4. Execute $WENOEXT/Allwmake to build the library 
