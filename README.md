![](documentation/WENOLogo.png)

![OpenFOAM Org](https://github.com/WENO-OF/WENOEXT/actions/workflows/c-ofORG.yml/badge.svg) 
![OpenFOAM ESI](https://github.com/WENO-OF/WENOEXT/actions/workflows/c-ofESI.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6589686.svg)](https://doi.org/10.5281/zenodo.6589686)



# WENO framework

Weighted essentially non-oscillatory library for the framework of OpenFOAM.
Detailed information about the theoretical background and the implementation can 
be found in:

 * [J. W. Gärtner, A. Kronenburg, T. Martin, Efficient WENO library for OpenFOAM, SoftwareX, 2020](./documentation/Gaertner2020.pdf)
 * [T. Martin and I. Shevchuk, Implementation and Validation of Semi-Implicit WENO Schemes Using OpenFOAM, Computation, 2018](./documentation/Martin2018.pdf)
 * [Solving the Level Set Equation using High-order Non-oscillatory Reconstruction](./documentation/Martin_Solving_the_Level_Set_Equation_using_High-order_Non-oscillatory_Reconstruction.pdf)

A quick overview of the WENO scheme is provided in this presentation:

 * [Presentation: WENOExt](./documentation/WENO-Presentation.pdf)

Please also check out the ***Discussion*** board of GitHub for more information about features or if you want to post a new idea.

**Versions:**

Major development stages of the library are marked by tags and recently also have a release with a DOI. 

### Supported OpenFOAM Versions:

 * OpenFOAM (ESI) v2406 (tested in CI on Ubuntu 24.04)

## Authors

 * Tobias Martin <tobias.martin@ntnu.no>
 * Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

When using this work please cite:
> J. W. Gärtner, A. Kronenburg, and T. Martin, “SoftwareX Efficient WENO library for OpenFOAM”, SoftwareX, vol. 12, p. 100611, 2020, doi: 10.1016/j.softx.2020.100611.

> T. Martin, and I. Shevchuk, “SoftwareX Efficient WENO library for OpenFOAM,” Computation, vol. 6(6), 2018, doi: 10.3390/computation6010006.

## Installation

1. Clone the directory with
    `git clone https://github.com/WENO-OF/WENOEXT.git`

2. Install the [Eigen](https://eigen.tuxfamily.org) library (version 3.4 or newer),
   e.g. on Ubuntu with `sudo apt install libeigen3-dev`

3. Execute `Allwmake` to build the library
   The compilation of the library uses cmake instead of wmake files!

### Options

Parallel compilation can be activated with:
`./Allwmake -j <# of cores>`
if the number of cores is omitted, it is determined automatically.

On the master branch only clean git commits can be compiled. 
If a dirty git state shall be compiled the `-f|--force` option has to be used
on the master branch. 

### CMake Options

Several options can be set over CMake commands with,
```bash
./Allwmake -D<CMAKE_OPTION>
```

|Cmake Command |Options| Description |
|:-----------------------|:---------|:----------------|
|MARCH_NATIVE    |ON/OFF | Activates `march=native` flag. Default ON </br> Use this flag if you get an "illegal instruction error" during execution.|
|USE_FMA|ON/OFF|Use std::fma for WENO math functions. Default ON|
|CMAKE_BUILD_TYPE|Release/Debug/None|When the debug option is selected the OpenFOAM FULLDEBUG flag is activated|

Commands not listed in the table are forwarded to cmake, allowing to use all standard CMake commands.

### Note to GNU compiler:

GNU compiler version must be higher than 7. For g++ < v7 an error is reported for 
the specialisation template syntax. 
The syntax in the code is according to C++11 standard which is available for g++ v7 and higher. 
 

## Usage

To use the WENO scheme you have to add the library to your controlDict by editing `system/controlDict`

    libs("libWENOEXT.so")

Within your `system/fvSchemes` file,

    divSchemes
    {
        div(phi,U) 	Gauss WENOUpwindFit 2 1;
    }

Here the first index '2' represents the order of the WENO scheme and the second index can be either
'1' for bounded or '0' for unbounded.

Further options can be set in the WENODict located in the 'system/' folder:

```C++
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      WENODict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// This dict contains expert parameters, which modify the standard WENO scheme.

    //- Stencil extension ratio:
    //  - < 2.5 :   decreased computational effort. May also decrease stability
    //              and accuracy
    //  - > 2.5 :   higher stability. May influence the accuracy of the SVD
    extendRatio     2.5;

    //- WENO stencil weighting parameters:
    p               4.0;
    dm              1000.0;

    //- Calculate best conditioned matrix
    //  This can save memory especially for high order WENO scheme
    //  Increases the calculation time! Default is off
    bestConditioned true;
    
    writeData       true; // Write out the collected stencil list and matrix data
                          // default is 'true' 

    maxCondition    1E-05;// Inverse of the maximum condition that the pseudo 
                          // inverse can have. Only change if you know what you
                          // are doing!

    checkCondition  false;// Check the condition of the pseudo inverse matrix
                          // If the central stencil has at least one zero entry
                          // the matrix is removed for all stencils of this cell.
// ************************************************************************* /
```

### Specialized Version for Scalar Transport

The limited WENOUpwindFit scheme uses a cell limited approach known from other
schemes such as linearUpwind. For scalar transport in the range 0 to 1 a 
specialized scheme called *WENOUpwindFit01* is available. This schemes limits
the value using the limiter of Zhang and Shu.




## Tutorials

The code contains two tutorials from the standard cavity test. 
To run the tutorials execute `./Allrun` in the tutorial/ directory.

## Tests

Testing is performed with the CATCH2 framework. You can compile and execute the tests
by executing `./runTest` in the test directory. Further instructions are found in [here](tests/TestInstructions.md) 

## Contineous Integration

To check the code, it is copied into a docker container created with the 
Dockerfiles in `CI/`. In the docker container the code is compiled
and unit tests executed. 

Executing the tests is controlled over the Makefile in the root directory:

```
# Just executing make runs all tests
make

# For OpenFOAM v2406
make runTestsOFv2406
```

This implementation allows to have the same testing on the local and remote branch.


## License 

This OpenFOAM library is under the GNU General Public License. 
The matrix operations use the header-only [Eigen](https://eigen.tuxfamily.org) library, 
which is licensed under the Mozilla Public License 2.0 (MPL2) and is not distributed with this library.



