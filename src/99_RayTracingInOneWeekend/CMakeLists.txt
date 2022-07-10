cmake_minimum_required(VERSION 3.16)
project(InWeekend)

message(${CMAKE_CURRENT_SOURCE_DIR})

# C++ 17
set(CMAKE_CXX_STANDARD 17)

# include
set(INCLUDE_H ${PROJECT_SOURCE_DIR}/include)
include_directories(${INCLUDE_H})

# list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# # Libigl
# option(LIBIGL_GLFW "Build target igl::glfw" ON)
# include(libigl)

# sub dir
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/01_create_image)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/02_vec3_class)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/03_ray_color)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/04_hit_sphere)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/05_hittable_objects)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/06_antialiasing)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/07_diffuse_materials)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/08_metal_materials)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/09_Dielectrics)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/10_positionable_camera)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/11_depth_of_field)


# # Add your project files
# file(GLOB SRC_FILES *.cpp)
# add_executable(${PROJECT_NAME} ${SRC_FILES})
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
