# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)

find_package(CGAL REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
set(LIBIGL_INCLUDE_DIR "C:/Users/Rached/Documents/Telecom/igd/test/Mesh Coloring/libigl/include")
list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../cmake")
include(libigl)
# libigl







#file(GLOB SRCFILES *.cpp)
#add_executable(my_executable ${SRCFILES})


# create a target per cppfile
add_executable(my_executable Source.cpp Operators.cpp HalfedgeDS.cpp HalfedgeBuilder.cpp)
target_link_libraries(my_executable CGAL::CGAL igl::core igl::opengl_glfw)
