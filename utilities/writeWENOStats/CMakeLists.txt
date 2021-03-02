# CMake File to Create the Library


add_executable(writeWENOStats 
    writeWENOStats.C
)




target_include_directories(writeWENOStats PUBLIC
    WENOEXT 
)

target_link_libraries(writeWENOStats PUBLIC
 WENOEXT
 -L$ENV{FOAM_LIBBIN}
)

set_target_properties(writeWENOStats PROPERTIES LINK_FLAGS "-fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

install(
    TARGETS writeWENOStats 
    DESTINATION $ENV{FOAM_USER_APPBIN} 
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE
)
