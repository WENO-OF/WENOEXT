![](documentation/WENOLogo.png)
![OpenFOAM 5.x](https://github.com/WENO-OF/WENOEXT/actions/workflows/c-of5.yml/badge.svg) 
![OpenFOAM 7](https://github.com/WENO-OF/WENOEXT/actions/workflows/c-of7.yml/badge.svg) 
![OpenFOAM 8](https://github.com/WENO-OF/WENOEXT/actions/workflows/c-of8.yml/badge.svg)

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

Different versions of the code are structured through tags:

 * [2.0] Improved speed by using Blaze matrix operations
 * [1.0] Version for OpenFOAM 5.x, 6 and 7.
             Uses improved parallising and C++11 features 
 * [0.1] Version for OpenFOAM 2.3.x to OpenFOAM 5.x 

### Supported OpenFOAM Versions:

 * OpenFOAM (ORG) v5.x - 8
 * OpenFOAM (ESI) v1912

## Authors

 * Tobias Martin <tobias.martin@ntnu.no>
 * Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

When using this work please cite:
> J. W. Gärtner, A. Kronenburg, and T. Martin, “SoftwareX Efficient WENO library for OpenFOAM”, SoftwareX, vol. 12, p. 100611, 2020, doi: 10.1016/j.softx.2020.100611.

> T. Martin, and I. Shevchuk, “SoftwareX Efficient WENO library for OpenFOAM,” Computation, vol. 6(6), 2018, doi: 10.3390/computation6010006.

## Installation

1. Clone the directory with
    `git clone https://github.com/WENO-OF/WENOEXT.git`

2. Execute `Allwmake` to build the library
   The compilation of the library uses cmake instead of wmake files!


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

## Tutorials

The code contains two tutorials from the standard cavity test. 
To run the tutorials execute `./Allrun` in the tutorial/ directory.

## Tests

Testing is performed with the CATCH2 framework. You can compile and execute the tests
by executing `./runTest` in the test directory. Further instructions are found in [here](tests/TestInstructions.md) 

## Contineous Integration

To check the code for different OpenFOAM implementations, e.g. OpenFOAM 8,
OpenFOAM v1912, etc. the code is copied into different docker containers
created with the Dockerfiles in `CI/`. In the docker container the code is compiled
and unit tests executed. 

Executing the tests is controlled over the Makefile in the root directory. The
different OpenFOAM versions can be tested with 

```
# Just executing make runs all tests
make

# For OpenFOAM v5.x
make runTestsOF5

# For OpenFOAM 8
make runTestsOF8
```

This implementation allows to have the same testing on the local and remote branch.


## License 

This OpenFOAM library is under the GNU General Public License. This library contains the [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/) library licensed under the BSD license. 
Redistribution and use of the Blaze source code with or without modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
  * Neither the names of the **Blaze** development group nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.



