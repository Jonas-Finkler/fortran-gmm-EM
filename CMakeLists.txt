cmake_minimum_required(VERSION 3.9)
project(heckdarting Fortran)
enable_language(Fortran)


option(INTEL "compile with ifort and mkl" OFF)
option(OPENMP "compile with openmp support" OFF)


if(CLUSTER)
    if(INTEL)
    else()
        set(CMAKE_Fortran_COMPILER /scicore/soft/apps/OpenMPI/3.1.2-GCC-8.2.0-2.31.1/bin/mpifort)
    endif()
else()
endif()

if(INTEL)
    set(CMAKE_Fortran_COMPILER ifort)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -ipo -mkl -Ofast ")
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -g -traceback -O2")
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
    endif()
else()
    set(CMAKE_Fortran_COMPILER gfortran)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -march=native -flto")
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -Wconversion -g -fbacktrace")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()

add_executable(
        em
        src/emExample.f90
        src/gmmEM.f90
        src/linalg.f90
        src/precision.f90
        src/random.f90
)

if(NOT INTEL)
    target_link_libraries(em ${libraries} openblas)
endif()

