project(MS_Intersection)

cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    set(OpenMP_C_LIB_NAMES "libomp")
    set(OpenMP_CXX_LIB_NAMES "libomp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/include")
    set(OpenMP_omp_LIBRARY "/usr/local/lib/libomp.dylib")
    set(OpenMP_libomp_LIBRARY "/usr/local/lib/libomp.dylib")
endif(APPLE)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(MS_Intersection
    src/surface/atom_sphere.h
    src/surface/atom_sphere.cpp
    src/surface/concave_edge.h
    src/surface/concave_edge.cpp
    src/surface/concave_sphere.h
    src/surface/concave_sphere.cpp
    src/surface/saddle_face.h
    src/surface/saddle_face.cpp
    src/surface/torus.h
    src/surface/torus.cpp

    src/grid/grid_edge.h
    src/grid/grid_edge.cpp
    src/grid/grid_point.h
    src/grid/grid_point.cpp
    
    src/utility/Jenkins_Traub.h
    src/utility/Jenkins_Traub.cpp
    src/utility/Matrix44.h
    src/utility/Matrix44.cpp
    src/utility/Types.h
    src/utility/Vector3d.h
    src/utility/Vector3d.cpp

    src/parser/pqr_parser.h
    src/parser/pqr_parser.cpp

    src/visualization/space_fill.h
    src/visualization/space_fill.cpp
    src/visualization/surface_patch.h
    src/visualization/surface_patch.cpp
    src/visualization/intersection_visual.h
    src/visualization/intersection_visual.cpp
    src/visualization/chain_visual.h
    src/visualization/chain_visual.cpp
    src/visualization/parallel_visual.h
    src/visualization/parallel_visual.cpp

    src/intersection/molecular_surface.h
    src/intersection/molecular_surface.cpp
    src/intersection/parallel_wrapper.h
    src/intersection/parallel_wrapper.cpp

    src/marchingcubes/LevelSet.h
    src/marchingcubes/LevelSet.cpp
    src/marchingcubes/LookUpTable.h
    src/marchingcubes/MarchingCubes.h
    src/marchingcubes/MarchingCubes.cpp

    src/main.cpp
    
    )
target_link_libraries(MS_Intersection PRIVATE OpenMP::OpenMP_CXX)