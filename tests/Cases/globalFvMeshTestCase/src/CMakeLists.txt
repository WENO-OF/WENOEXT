# CMake File to Create the Library



add_executable(globalFvMeshTest.exe 
    globalFvMesh-Test.C
)




target_include_directories(globalFvMeshTest.exe PUBLIC
    WENOEXT
    $ENV{FOAM_SRC}/sampling/lnInclude
)

target_link_libraries(globalFvMeshTest.exe PUBLIC
 WENOEXT
 -L$ENV{FOAM_LIBBIN}
 -lsampling
 -ldl
 -lm 
)

set_target_properties(globalFvMeshTest.exe PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

install(
    TARGETS globalFvMeshTest.exe 
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
