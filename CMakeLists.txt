# ===========================================================================
#       ██╗    ██╗███████╗███╗   ██╗ ██████╗     ███████╗██╗  ██╗████████╗
#       ██║    ██║██╔════╝████╗  ██║██╔═══██╗    ██╔════╝╚██╗██╔╝╚══██╔══╝
#       ██║ █╗ ██║█████╗  ██╔██╗ ██║██║   ██║    █████╗   ╚███╔╝    ██║   
#       ██║███╗██║██╔══╝  ██║╚██╗██║██║   ██║    ██╔══╝   ██╔██╗    ██║   
#       ╚███╔███╔╝███████╗██║ ╚████║╚██████╔╝    ███████╗██╔╝ ██╗   ██║   
#        ╚══╝╚══╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝     ╚══════╝╚═╝  ╚═╝   ╚═╝   
# ===========================================================================                                                                                                                                                        
# License
#    This file is part of WENO Ext.
#
#    WENO Ext is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    WENO Ext is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with  WENO Ext.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Author
#   Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020
#
# ===========================================================================
# CMake File for WENOEXT Project
# CMake is used instead of OpenFOAM standard wmake as it allows 
# for some more flexibility between the different OpenFOAM distributions
# without writing extensive bash scripts to do this 
# Further, it integrates directly with Catch2 and the Blaze library if
# needed
# ===========================================================================

cmake_minimum_required (VERSION 3.10)
# make cache variables for install destinations
include(GNUInstallDirs)


# Enforce C++ 14 standard
# C++ 14 is required for blaze
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Set minimum gnu compiler version
set(GNU_MIN_VERSION 7.5.0)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (NOT CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL ${GNU_MIN_VERSION})
        message(WARNING "GNU compiler version ${CMAKE_CXX_COMPILER_VERSION} is below ${GNU_MIN_VERSION}")
    endif()
endif()



# Check that OpenFOam environment variables are available
if (NOT DEFINED ENV{FOAM_USER_LIBBIN} OR NOT DEFINED ENV{FOAM_SRC})
    message(FATAL_ERROR "OpenFOAM environemt variables not available")
endif()

#===============================================================================
#   Switch to activate march=native 
#===============================================================================
# If the code is distributed on several different CPU architectures swtich
# this option off
# Default is ON

option (MARCH_NATIVE "If enabled, uses -march=native" ON)

# Blaze requires LAPACK for the calculation of SVD and the eigenvalues of 
# a matrix.
# If Blaze functions should be used it has to be switched on
# Default: OFF
option (USE_LAPACK "If enabled, uses LAPACK for matrix operations" OFF)

# If AVX support is enabled set the FP_FAST_FMA flag to use
# std::fma instead of (a*b) + c for the mathFunctionsWENO::det2() function 
option (USE_FMA "If enabled, uses fused multiply and add std::FMA" ON)

#===============================================================================
#   Create Version File
#===============================================================================

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

execute_process(
    COMMAND git describe --tags  
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Process tag string to set the CMake project version
string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" VERSION_MAJOR "${GIT_TAG}")
string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${GIT_TAG}")
string(REGEX REPLACE "^[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${GIT_TAG}")

if ("${VERSION_MINOR}" STREQUAL "")
    set(VERSION_MINOR "0")
endif()
if ("${VERSION_MAJOR}" STREQUAL "")
    set(VERSION_MAJOR "0")
endif()

set(WENOEXT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR})

configure_file(
    ${CMAKE_SOURCE_DIR}/versionRules/version.h.in
    ${CMAKE_BINARY_DIR}/generated/version.h
)

include_directories(${CMAKE_BINARY_DIR}/generated)


# Project name 
project (WENOEXT VERSION ${WENOEXT_VERSION})

# Set default install path
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/installed" CACHE PATH "default install path" FORCE )
endif()



#===============================================================================
#   Set OpenFOAM Environment
#===============================================================================


message("\n-------------------------------------------------------------------") 
message("\tGenerate WENO code rules for OpenFOAM version $ENV{WM_PROJECT_VERSION}")
message("-------------------------------------------------------------------")

# Check the OpenFOAM version and set the correct include statements
string(REGEX MATCH "([0-9]+)" _ $ENV{WM_PROJECT_VERSION})
set(OF_VERSION_NUMBER ${CMAKE_MATCH_1})

# If VERSION > 1000 then it is a OpenFOAM.com (ESI) version
# In OpenFOAM (ESI) the triSurface is moved into the surfaceMesh 
# Note: The if() statement evaluates a set variable to true if no
# $-Sign is used. For the cmakedefine statement this means that all set 
# variables are resolved as true
if (${OF_VERSION_NUMBER} GREATER_EQUAL 1000)
    SET(OF_FORK_ESI true BOOL "OpenFOAM Fork: ESI")
    SET(OF_FORK_VERSION ${OF_VERSION_NUMBER})
    message("OF_VERSION: ${OF_FORK_VERSION}")
else()
    SET(OF_FORK_ORG true BOOL "OpenFOAM Fork: Org")
    SET(OF_FORK_VERSION ${OF_VERSION_NUMBER})
    message("OF_VERSION: ${OF_FORK_VERSION}")
endif()



configure_file(
    ${CMAKE_SOURCE_DIR}/versionRules/codeRules.H.in
    ${CMAKE_BINARY_DIR}/generated/codeRules.H
)

# Install generated files
install(FILES ${CMAKE_BINARY_DIR}/generated/version.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(FILES ${CMAKE_BINARY_DIR}/generated/codeRules.H
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")


# If VERSION > 1000 then it is a OpenFOAM.com (ESI) version
# In OpenFOAM (ESI) the triSurface is moved into the surfaceMesh 
if (${OF_VERSION_NUMBER} GREATER_EQUAL 1000)
    SET(TRISURFACE_INC)
    SET(TRISURFACE_LIB)
else()
    SET(TRISURFACE_INC "$ENV{FOAM_SRC}/triSurface/lnInclude")
    SET(TRISURFACE_LIB "-ltriSurface")
endif()




# Set c++ flags
add_compile_options(
    -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-deprecated
    -ftemplate-depth-100
    -O3
    )

if (${MARCH_NATIVE})
    add_compile_options(
        -march=native
    )
endif()

if (${USE_FMA})
    add_definitions( -DUSE_FMA )
endif()

set(CMAKE_SHARED_LINKER_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")


add_definitions(-Dlinux64 -DWM_ARCH_OPTION=$ENV{WM_ARCH_OPTION} -DWM_DP -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE} -DNoRepository)



#===============================================================================
#                           Add subdirectories
#===============================================================================

# Directory to build the library
add_subdirectory(blaze-3.8)
add_subdirectory(libWENOEXT)
add_subdirectory(utilities)



message("\n-------------------------------------------------------------------") 
message("\tBuild State")
message("-------------------------------------------------------------------")
if (CMAKE_BUILD_TYPE)
    message("Built state is: ${CMAKE_BUILD_TYPE}")
else()
    message("Built state is: unspecified")
endif()
# Build only tests if activated
# As strequal is a case sensitive comparison and the cmake build type can be 
# specified case insensitive it has to be first converted
string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
if (build_type_upper STREQUAL "RELEASE")
    message("Tests are NOT build for the release type")
else()
    message("Building tests")
    message("This can be deactivated by using CMAKE_BUILD_TYPE=Release")
    add_subdirectory(tests)
endif()
message("-------------------------------------------------------------------\n")


#===============================================================================
#                      Install CMake Exports
#===============================================================================

install(EXPORT WENOEXTTargets 
        FILE WENOEXT.cmake
        NAMESPACE WENOEXT::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/WENOEXT
)

