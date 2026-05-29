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

add_library(WENOEXT SHARED 
    BlazeIO/BlazeIO.C
    WENOBase/geometryWENO/geometryWENO.C
    WENOBase/WENOBase.C
    WENOBase/globalfvMesh.C
    WENOBase/matrixDB.C
    WENOBase/reconstructRegionalMesh.C
    WENOUpwindFit/makeWENOUpwindFit.C
    WENOUpwindFit01/makeWENOUpwindFit01.C
)



target_include_directories(WENOEXT PUBLIC
  $ENV{FOAM_SRC}/dynamicMesh/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/surfMesh/lnInclude
  $ENV{FOAM_SRC}/fileFormats/lnInclude 
  ${TRISURFACE_INC}
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/blaze-3.8>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/WENOBase/geometryWENO>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/WENOBase>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/BlazeIO>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Blaze>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libWENOEXT/WENOBase/geometryWENO>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libWENOEXT/WENOBase>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libWENOEXT/BlazeIO>"
)


# Set OpenFOAM library directory
set(FOAM_LIBBIN $ENV{FOAM_LIBBIN})

# Find all required OpenFOAM libraries
find_library(LIB_DYNAMIC_MESH dynamicMesh PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_FINITE_VOLUME finiteVolume PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_MESH_TOOLS meshTools PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_SPECIE specie PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_SURF_MESH surfMesh PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_FILE_FORMATS fileFormats PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_OPENFOAM OpenFOAM PATHS ${FOAM_LIBBIN} NO_DEFAULT_PATH)
find_library(LIB_PSTREAM Pstream PATHS ${FOAM_LIBBIN}/dummy NO_DEFAULT_PATH)

# Check if libraries were found
if(NOT LIB_OPENFOAM)
    message(FATAL_ERROR "OpenFOAM libraries not found in ${FOAM_LIBBIN}")
endif()

# Link libraries
target_link_libraries(WENOEXT PUBLIC
    Blaze
    ${LIB_DYNAMIC_MESH}
    ${LIB_FINITE_VOLUME}
    ${LIB_MESH_TOOLS}
    ${LIB_SPECIE}
    ${LIB_SURF_MESH}
    ${LIB_FILE_FORMATS}
    ${LIB_OPENFOAM}
    ${LIB_PSTREAM}
    ${TRISURFACE_LIB}
)

target_compile_definitions(WENOEXT PRIVATE
    $<$<CONFIG:Debug>:
        FULLDEBUG
    >
)



install(
    TARGETS WENOEXT 
    EXPORT WENOEXTTargets
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the header files in the include directory 

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.H"
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/WENOBase/WENOCoeff.C
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libWENOEXT/WENOBase/
)
