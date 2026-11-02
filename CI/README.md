# Dockerfiles for Contineous Integration

Dockerfiles to generate the docker image to compile the library and execute 
the unit tests with,

 * OpenFOAM v2406 (ESI) on Ubuntu 24.04 with Eigen 3.4

The image is self-contained: OpenFOAM is installed from the precompiled
openfoam.com Debian packages and Eigen from the Ubuntu package `libeigen3-dev`.

Build the image and run the tests with `make runTestsOFv2406` from the 
root directory.

For other ESI versions adapt the package name `openfoam2406-dev` and the 
`foamDotFile` path in `Dockerfile.OFv2406`.
