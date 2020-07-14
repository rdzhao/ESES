project(MS_Intersection)

cmake_minimum_required(VERSION 3.12)

set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/9.0.0_1/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/9.0.0_1/bin/clang++")
set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/9.0.0_1/lib")
set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/9.0.0_1/include")

OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

# Find OpenMP
if(APPLE AND USE_OpenMP)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

if(USE_OpenMP)
  find_package(OpenMP REQUIRED)
endif(USE_OpenMP)

if (OPENMP_FOUND)
    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

#aux_source_directory(. srcs)

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