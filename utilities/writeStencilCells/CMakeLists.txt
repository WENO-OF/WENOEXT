# CMake File to Create the Library


add_executable(writeStencilCells
    writeStencilCells.C
)




target_include_directories(writeStencilCells PUBLIC
    WENOEXT 
)

target_link_libraries(writeStencilCells PUBLIC
 WENOEXT
 -L$ENV{FOAM_LIBBIN}
)

set_target_properties(writeStencilCells PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

install(
    TARGETS writeStencilCells 
    DESTINATION $ENV{FOAM_USER_APPBIN} 
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
