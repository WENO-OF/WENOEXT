# CMake File to Create the Library


add_library(WENOEXT SHARED 
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


install(
    TARGETS WENOEXT 
    DESTINATION $ENV{FOAM_USER_LIBBIN}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
