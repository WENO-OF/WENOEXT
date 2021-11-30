#===============================================================================
#                           Blaze 3.8 CMake File
#===============================================================================
# The Blaze library includes a logging function that interferes with 
# OpenFOAMs own log function. Therefore, it is removed from the blaze library 
# here.
# Providing own CMakeLists file 
#===============================================================================

add_library(Blaze INTERFACE)

target_include_directories(Blaze INTERFACE  .)

#===============================================================================
# Configure LAPACK
#===============================================================================

# Blaze requires LAPACK for the calculation of SVD and the eigenvalues of 
# a matrix.

find_package(LAPACK REQUIRED)
target_link_libraries(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LIBRARIES}>)
target_compile_options(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LINKER_FLAGS}>)


#===============================================================================
#   Set Cache of Blaze Automatic 
#===============================================================================

message("\n-------------------------------------------------------------------") 
message("\tSet blaze cache size")
message("-------------------------------------------------------------------")
message("Cache size is automatically set if not provided by the user with")
message("\tcmake ../ -DBLAZE_CACHE_SIZE=<Cache Size>")
if(NOT BLAZE_CACHE_SIZE)
    set(BLAZE_CACHE_SIZE "auto")
endif()


    if ("${BLAZE_CACHE_SIZE}" STREQUAL "auto")
        message("Automatic Cache Size Configuration")
        set(flag 1)

        execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index3/size
                   OUTPUT_VARIABLE tmp
                   RESULT_VARIABLE flag
                   ERROR_QUIET)
        if (flag)
            execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index2/size
                     OUTPUT_VARIABLE tmp
                     RESULT_VARIABLE flag
                     ERROR_QUIET)
        endif (flag)
        if (flag)
            execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index1/size
                     OUTPUT_VARIABLE tmp
                     RESULT_VARIABLE flag
                     ERROR_QUIET)
        endif (flag)

         

        if (flag)
           message(WARNING "Cache size not found automatically. Using default value as cache size.")
           set(tmp "3072")
        endif (flag)

        string(REGEX MATCH "([0-9][0-9]+)" tmp ${tmp}) # Get a number containing at least 2 digits in the string tmp
        math(EXPR BLAZE_CACHE_SIZE ${tmp}*1024) # Convert to bytes (assuming that the value is given in kibibytes)
        set(BLAZE_CACHE_SIZE "${BLAZE_CACHE_SIZE}UL")
        target_compile_definitions( Blaze INTERFACE BLAZE_CACHE_SIZE=${BLAZE_CACHE_SIZE} )
        message("Configuring Cache Size : ${BLAZE_CACHE_SIZE}")
    elseif(${BLAZE_CACHE_SIZE})
        message("Blaze Cache Size manually set")
        target_compile_definitions( Blaze INTERFACE BLAZE_CACHE_SIZE=${Blaze_CACHE_SIZE} )
        message("Configuring Cache Size : ${BLAZE_CACHE_SIZE}")
    else()
        message("Using default Cache Size.")
    endif ()

message("-------------------------------------------------------------------")
unset(BLAZE_CACHE_SIZE CACHE)
