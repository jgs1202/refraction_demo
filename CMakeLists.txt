cmake_minimum_required(VERSION 3.13)
project(refraction)

set(CMAKE_CXX_STANDARD 14)

add_executable(refraction main.cpp)

INCLUDE_DIRECTORIES (/System/Library/Frameworks)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
target_link_libraries(refraction "/System/Library/Frameworks/GLUT.framework" "/System/Library/Frameworks/openGL.framework")