#===============================================================================
#                           Blaze 3.8 CMake File
#===============================================================================
# The Blaze library includes a logging function that interferes with 
# OpenFOAMs own log function. Therefore, it is removed from the blaze library 
# here.
# Providing own CMakeLists file 
#===============================================================================

add_library(Blaze INTERFACE)

target_include_directories(Blaze INTERFACE
   "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
   "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Blaze>"
)

install(TARGETS Blaze
        EXPORT WENOEXTTargets 
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Blaze)

install(DIRECTORY blaze DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Blaze)




#===============================================================================
# Configure LAPACK
#===============================================================================

if (USE_LAPACK)
    find_package(LAPACK REQUIRED)
    target_link_libraries(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LIBRARIES}>)
    target_compile_options(Blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LINKER_FLAGS}>)
endif()


#===============================================================================
#   Set Cache of Blaze Automatic 
#===============================================================================

message("\n-------------------------------------------------------------------") 
message("\tSet blaze cache size")
message("-------------------------------------------------------------------")
message("Cache size is automatically set if not provided by the user with")
message("\tcmake ../ -DBLAZE_CACHE_SIZE=<Cache Size>")

set(_blaze_default_cache_kib 3072)

if(NOT BLAZE_CACHE_SIZE OR BLAZE_CACHE_SIZE STREQUAL "auto")
  message("Automatic Cache Size Configuration")
  set(_best_level 0)
  set(_best_kib "")

  # Pick the highest-level data/unified cache of cpu0 by reading level/type/size,
  # instead of assuming index3 == L3 (index1 is usually the L1 *instruction* cache).
  file(GLOB _cache_dirs "/sys/devices/system/cpu/cpu0/cache/index*")
  foreach(_dir IN LISTS _cache_dirs)
    if(EXISTS "${_dir}/level" AND EXISTS "${_dir}/type" AND EXISTS "${_dir}/size")
      file(READ "${_dir}/level" _lvl)
      file(READ "${_dir}/type"  _type)
      file(READ "${_dir}/size"  _size)
      string(STRIP "${_lvl}"  _lvl)
      string(STRIP "${_type}" _type)
      string(STRIP "${_size}" _size)
      if(_type STREQUAL "Instruction")
        continue()
      endif()
      if(_size MATCHES "^([0-9]+)([KMG]?)$")
        set(_kib ${CMAKE_MATCH_1})
        if(CMAKE_MATCH_2 STREQUAL "M")
          math(EXPR _kib "${_kib} * 1024")
        elseif(CMAKE_MATCH_2 STREQUAL "G")
          math(EXPR _kib "${_kib} * 1024 * 1024")
        elseif(CMAKE_MATCH_2 STREQUAL "")
          math(EXPR _kib "${_kib} / 1024")
        endif()
        if(_lvl GREATER _best_level AND _kib GREATER 0)
          set(_best_level ${_lvl})
          set(_best_kib ${_kib})
        endif()
      endif()
    endif()
  endforeach()

  # Fallback: glibc's getconf (returns bytes, may be 0/empty)
  if(NOT _best_kib)
    foreach(_lv 3 2 1)
      execute_process(COMMAND getconf LEVEL${_lv}_CACHE_SIZE
                      OUTPUT_VARIABLE _bytes OUTPUT_STRIP_TRAILING_WHITESPACE
                      RESULT_VARIABLE _rc ERROR_QUIET)
      if(_rc EQUAL 0 AND _bytes MATCHES "^[0-9]+$" AND _bytes GREATER 0)
        math(EXPR _best_kib "${_bytes} / 1024")
        set(_best_level ${_lv})
        break()
      endif()
    endforeach()
  endif()

  if(NOT _best_kib)
    message(WARNING "Cache size not found automatically. Using default value as cache size.")
    set(_best_kib ${_blaze_default_cache_kib})
  else()
    message("Detected L${_best_level} cache: ${_best_kib} KiB")
  endif()

  math(EXPR _bytes "${_best_kib} * 1024")
  # Blaze (system/CacheSize.h) asserts 100000 < cacheSize < 100000000 bytes
  if(_bytes GREATER_EQUAL 100000000)
    message("Detected cache exceeds Blaze's upper limit, clamping to 99999744 bytes")
    set(_bytes 99999744)   # largest multiple of 1024 below 1e8
  elseif(_bytes LESS_EQUAL 100000)
    message("Detected cache below Blaze's lower limit, clamping to 102400 bytes")
    set(_bytes 102400)
  endif()
  target_compile_definitions(Blaze INTERFACE BLAZE_CACHE_SIZE=${_bytes}UL)
  message("Configuring Cache Size : ${_bytes}UL")

elseif(BLAZE_CACHE_SIZE MATCHES "^([0-9]+)(UL)?$")
  message("Blaze Cache Size manually set")
  target_compile_definitions(Blaze INTERFACE BLAZE_CACHE_SIZE=${CMAKE_MATCH_1}UL)
  message("Configuring Cache Size : ${CMAKE_MATCH_1}UL")
else()
  message(WARNING "Invalid BLAZE_CACHE_SIZE='${BLAZE_CACHE_SIZE}' (expected bytes). Using Blaze's built-in default.")
endif()

message("-------------------------------------------------------------------")
unset(BLAZE_CACHE_SIZE CACHE)
