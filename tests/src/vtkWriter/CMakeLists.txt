# CMake File to create vtkWriter library
cmake_minimum_required(VERSION 3.10)


file(GLOB Sources "*.C")

add_library(vtkWriter STATIC ${Sources})

# Required for OpenFOAM
target_include_directories(vtkWriter PUBLIC
    $ENV{FOAM_SRC}/finiteVolume/lnInclude
    $ENV{FOAM_SRC}/meshTools/lnInclude
    $ENV{FOAM_SRC}/OpenFOAM/lnInclude
    $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
)

target_link_libraries(vtkWriter PUBLIC
    -lfiniteVolume
    -lOpenFOAM
    -lmeshTools
)
