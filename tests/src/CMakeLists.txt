# CMake File to Create the Library


add_subdirectory(KDTree)

add_subdirectory(vtkWriter)

add_executable(WENO_TEST 
    main.C
    geometryWENO-BasicFunc-Test.C
    mathFunctionsWENO-Test.C
    matrixDB-Test.C
    WENOUpwindFit-Test.C
    WENOUpwindFit-Test-parallel.C
    WENOUpwindFit-AdvectionTest.C 
    WENOBaseIO-Test.C
    List3D-Test.C
)




target_include_directories(WENO_TEST PUBLIC
    WENOEXT
    KDTree
    vtkWriter
    ./Catch2SingleHeader
    $ENV{FOAM_SRC}/sampling/lnInclude
)

target_link_libraries(WENO_TEST PUBLIC
 WENOEXT
 KDTree
 vtkWriter
 -L$ENV{FOAM_LIBBIN}
 -lsampling
 -ldl
 -lm 
)

set_target_properties(WENO_TEST PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

install(
    TARGETS WENO_TEST 
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
