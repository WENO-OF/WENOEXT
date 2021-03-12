# CMake File to Create the Library



add_executable(performanceRun.exe 
    main.C
)




target_include_directories(performanceRun.exe PUBLIC
    WENOEXT
    $ENV{FOAM_SRC}/sampling/lnInclude
)

target_link_libraries(performanceRun.exe PUBLIC
 WENOEXT
 -L$ENV{FOAM_LIBBIN}
 -lsampling
 -ldl
 -lm 
)

set_target_properties(performanceRun.exe PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

install(
    TARGETS performanceRun.exe 
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
