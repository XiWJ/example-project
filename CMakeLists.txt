cmake_minimum_required(VERSION 3.16)
project(example)

# C++ 17
set(CMAKE_CXX_STANDARD 17)

# include
set(INCLUDE_H ${PROJECT_SOURCE_DIR}/include)
include_directories(${INCLUDE_H})

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

# sub dir
add_subdirectory(${PROJECT_SOURCE_DIR}/src/01_Normals)

# # Add your project files
# file(GLOB SRC_FILES *.cpp)
# add_executable(${PROJECT_NAME} ${SRC_FILES})
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
