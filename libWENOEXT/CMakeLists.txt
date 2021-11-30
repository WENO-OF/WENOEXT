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
  Blaze
  $ENV{FOAM_SRC}/dynamicMesh/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/surfMesh/lnInclude
  $ENV{FOAM_SRC}/fileFormats/lnInclude 
  ${TRISURFACE_INC}
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  ./WENOBase/geometryWENO
  ./WENOBase/
  ./BlazeIO/
  ../versionRules
)

target_link_libraries(WENOEXT PUBLIC
 Blaze
 -L$ENV{FOAM_LIBBIN}
 -ldynamicMesh 
 -lfiniteVolume 
 -lmeshTools
 -lspecie
 -lsurfMesh
 -lfileFormats
 -lOpenFOAM
 ${TRISURFACE_LIB}
)


target_compile_definitions(WENOEXT PRIVATE
    $<$<CONFIG:Debug>:
        FULLDEBUG
    >
)



install(
    TARGETS WENOEXT 
    DESTINATION $ENV{FOAM_USER_LIBBIN}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
