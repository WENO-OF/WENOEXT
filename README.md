WENO framework
====================
Weighted essentially non-oscillatory library for the framework of OpenFOAM.

Currently supported version:
	- 2.3.x
	

Installation
============
1. rename the folder in "WENOEXT"
2. move the folder "WENOEXT" into $FOAM_SRC
3. add the following line to your ~/.bashrc after ". ~/OpenFOAM/OpenFOAM-2.3.x/bashrc":
     . $FOAM_SRC/WENOEXT/bashrc
4. parse your ~/.bashrc or open a new terminal
5. Execute $WENOEXT/Allwmake to build the library 
