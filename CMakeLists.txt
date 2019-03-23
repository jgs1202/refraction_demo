cmake_minimum_required(VERSION 3.13)
project(refraction)

set(CMAKE_CXX_STANDARD 14)

add_executable(refraction main.cpp)

INCLUDE_DIRECTORIES (/System/Library/Frameworks)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
link_directories(${OpenCL_LIBRARY})

target_include_directories (opencl_sample PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(refraction "/System/Library/Frameworks/GLUT.framework" "/System/Library/Frameworks/openGL.framework" ${OpenCL_LIBRARY})