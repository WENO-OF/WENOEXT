# WENO Test Cases

Test cases for individual functions (unit tests) as well as 
tests for the complete scheme (integration tests). All tests are 
managed with the software Catch2 and are compiled into one executable
the different tests can then be run individually by selecting the 
appropiate tags. To list all available tags use:

    tests.exe -t 

For further information about command line arguments in Catch2 use 

    tests.exe -?

## Included Tests

Currently these tests are available

### 1. geometryWENO Class

A test function for the Jacobian and the gaussQuad function is provided.
Requires a 3D mesh given in 'globalFvMeshTest/case-3D'

Run with `tests.exe [3D]` in globalFvMeshTest/case-3D directory.


### 2. WENOUpwindFit

An integration test to compare the result of the WENO scheme to the linear scheme
for the calculation of the divergence of a scalar field. 

Run with `tests.exe [2D]` in the Case directory.

### 3. GlobalFvMesh test case

Test if the mapping of global to local cellID and reverse is correct. As this test
is done in parallel it is not included in the Catch2 environment but uses 
FatalError statements to print out error messages

## Mesh Study

To generate a small mesh study of the implemented WENO scheme the script
`WENOEXT/tests/Case/runMeshStudy.sh` can be executed. It runs on the 2D mesh 
the WENOUpwindFit test case for a 50x50, 100x100, 150x150 mesh storing the results
in a file. The results can be plotted with the gnuplot script provided in the same 
directory. 


