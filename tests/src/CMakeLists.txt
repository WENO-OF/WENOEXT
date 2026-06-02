# CMake File to Create the Library


add_subdirectory(KDTree)

add_subdirectory(vtkWriter)

add_executable(WENO_TEST 
    main.C
    geometryWENO-BasicFunc-Test.C
    mathFunctionsWENO-Test.C
    matrixDB-Test.C
    WENOUpwindFit-Test.C
    WENOUpwindFit-AdvectionTest.C 
    WENOBaseIO-Test.C
    List3D-Test.C
    globalFvMesh-Test.C
)




target_include_directories(WENO_TEST PUBLIC
    WENOEXT
    KDTree
    vtkWriter
    $ENV{FOAM_SRC}/sampling/lnInclude
)

target_link_libraries(WENO_TEST PUBLIC
 Catch2::Catch2
 WENOEXT
 KDTree
 vtkWriter
 -L$ENV{FOAM_LIBBIN}
 -lsampling
 -ldl
 -lm 
)


if (APPLE)
    set_target_properties(WENO_TEST PROPERTIES
        MACOSX_RPATH ON
        BUILD_RPATH "${CMAKE_BINARY_DIR}/lib;${CMAKE_BINARY_DIR}"
        INSTALL_RPATH "@executable_path/../lib;@loader_path/../lib"
    )
else()
set_target_properties(WENO_TEST PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")
endif()

install(
    TARGETS WENO_TEST 
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
