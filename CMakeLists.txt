project(MS_Intersection)

cmake_minimum_required(VERSION 2.8)

find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

file(GLOB Headers "*.h")
file(GLOB Sources "*.cpp")

add_executable(MS_Intersection ${Headers} ${Sources})