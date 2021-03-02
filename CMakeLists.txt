# CMake File for WENOEXT Project
# CMake is used instead of OpenFOAM standard wmake as it allows 
# for some more flexibility between the different OpenFOAM distributions
# without writing extensive bash scripts to do this 
# Further, it integrates directly with Catch2 and the Blaze library if
# needed


cmake_minimum_required (VERSION 3.10)

# Project name 
project (WENOEXT)

# Enforce C++ 14 standard
# C++ 14 is required for blaze
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Check that OpenFOam environment variables are available
if (NOT DEFINED ENV{FOAM_USER_LIBBIN} OR NOT DEFINED ENV{FOAM_SRC})
    message(FATAL_ERROR "OpenFOAM envrionemt variables not available")
endif()

# Get current git commit
execute_process(
  COMMAND  git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
  ${CMAKE_SOURCE_DIR}/versionRules/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
)

include_directories(${CMAKE_BINARY_DIR}/generated)



# Check the OpenFOAM version and set the correct include statements
string(REGEX MATCH "([0-9]*)" _ $ENV{WM_PROJECT_VERSION})
set(OF_VERSION_NUMBER ${CMAKE_MATCH_1})

# If VERSION > 1812 then it is a OpenFOAM.com (ESI) version
# In OpenFOAM (ESI) the triSurface is moved into the surfaceMesh 
if (${OF_VERSION_NUMBER} GREATER_EQUAL 1812)
    SET(TRISURFACE_INC)
    SET(TRISURFACE_LIB)
else()
    SET(TRISURFACE_INC "$ENV{FOAM_SRC}/triSurface/lnInclude")
    SET(TRISURFACE_LIB "-ltriSurface")
endif()

# Generate OpenFOAM Code rules 
message("Generate WENO code rules for OpenFOAM version $ENV{WM_PROJECT_VERSION}")
execute_process(COMMAND ${CMAKE_SOURCE_DIR}/versionRules/genCodeRules.sh ${CMAKE_SOURCE_DIR})


# Set c++ flags
set(CMAKE_CXX_FLAGS "-m64 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3 -ftemplate-depth-100 -Wno-deprecated")

set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

# Get the label size
#execute_process(COMMAND grep -rin "export WM_LABEL_SIZE" $ENV{WM_PROJECT_DIR}/etc COMMAND  sed  "s@.*export WM_LABEL_SIZE=@@" COMMAND  tr -d '\n' OUTPUT_VARIABLE LABEL_SIZE)
add_definitions(-Dlinux64 -DWM_ARCH_OPTION=$ENV{WM_ARCH_OPTION} -DWM_DP -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE} -DNoRepository)


# Directory to build the library
add_subdirectory(libWENOEXT)
add_subdirectory(utilities)
add_subdirectory(tests/src)

# Add here later for testing 



