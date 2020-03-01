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
Requires a 3D mesh given in 'case-3D'

### 2. WENOUpwindFit

An integration test to compare the result of the WENO scheme to the linear
and upwind scheme for a steady state convection diffusion problem. A comparison
to an analytical solution is given.

The test checks that the WENO scheme does not perform worse than the upwind scheme
and prints at the end the mean error to the analytical solution. 

** Note: For the WENO scheme at least 2 iterations have to be performed as the WENO
         explicit correction needs the updated field information. **

### 3. GlobalFvMesh test case

Test if the mapping of global to local cellID and reverse is correct. As this test
is done in parallel it is not included in the Catch2 environment but uses 
FatalError statements to print out error messages
