cmake_minimum_required(VERSION 2.8)
#project(KDTree)

#set(CMAKE_BUILD_TYPE Debug)

include_directories(.)
add_library(KDTree STATIC "KDTree.cpp")

#add_executable(test_kdtree "test_kdtree.cpp")
#target_link_libraries(test_kdtree KDTree)
