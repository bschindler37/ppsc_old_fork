
###################################################################################
#
# PPSC: Pseudo Particle Strong Coupling
#
# Copyright (C) 2017 by H. U.R. Strand
#
###################################################################################

cmake_minimum_required(VERSION 3.0)
project("ppsc")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(COMMON_LIBS "")

# ~~ Add HDF5 ~~	## copied from nessi CMakeLists.txt
if (hdf5)
  #find_library(HDF5_LIB hdf5)
  find_package(HDF5)
  set(HDF5_LIB ${LIBS} ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIR})
  #message(STATUS "${HDF5_INCLUDE_DIR}")
  #list(APPEND COMMON_LIBS ${HDF5_LIB})
  #message(STATUS "${HDF5_LIB}")
  add_definitions("-DCNTR_USE_HDF5")
endif (hdf5)

# ~~ OpenMP and MPI flags ~~

option(omp "Build with OpenMP support" OFF)
if (omp)
    message(STATUS "Building with OpenMP")
    find_package(OpenMP REQUIRED)
    add_definitions("-DCNTR_USE_OMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

option(mpi "Build with OpenMPI support" OFF)

if (mpi)
    message(STATUS "Building with OpenMPI")
    find_package(MPI REQUIRED)
    add_definitions("-DCNTR_USE_MPI")
    set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_COMPILE_FLAGS
        "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS
        "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND COMMON_LIBS ${MPI_LIBRARIES})
endif ()

# ~~ Add BOOST ~~
find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#list(APPEND COMMON_LIBS ${Boost_LIBRARIES})

find_package(Eigen3 REQUIRED NO_MODULE)
add_definitions("-DCNTR_USE_HDF5")
find_library(CNTR_LIBRARY REQUIRED NAMES cntr)
find_library(CNTR_HDF5_LIBRARY REQUIRED NAMES cntr_hdf5)
set(CNTR_INCLUDE_DIR "/usr/local/include/")
set(CNTR_HDF5_INCLUDE_DIR "/usral/include/cntr/hdf5")
message("CNTR_LIB ${CNTR_LIBRARY} INCLUDE ${CNTR_INCLUDE_DIR}")
message("CNTR_H5_LIB ${CNTR_HDF5_LIBRARY} INCLUDE ${CNTR_HDF5_INCLUDE_DIR}")
find_package(OpenMP REQUIRED)

# ~~ Relocatable code for shared library ~~
set(GCC_RELOCATABLE_FLAGS "-fPIC")
add_definitions(${GCC_RELOCATABLE_FLAGS})

# ~~ Warning Flags ~~
# Re-enable whenever you feel strong enough.
# set(GCC_WARNING_FLAGS "-Wall -Wextra -pedantic")
# add_definitions(${GCC_WARNING_FLAGS})

# ~~ CMake Configuration ~~
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})
link_directories($)

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(ppsc)
add_subdirectory(programs)

# ~~ Include Guard ~~
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

