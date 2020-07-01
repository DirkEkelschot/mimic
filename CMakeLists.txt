# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)

# set the project name and version
project(ADAPT VERSION 1.0)


# specify the C++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

include(platforms/MacOSX_Dirk_home.cmake)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wno-dev)

#set(LIBRARY_OUTPUT_PATH lib/)

set(SRC adapt_output.cpp adapt_compute.cpp adapt_part_func.cpp adapt_schedule.cpp adapt_operations.cpp main.cpp)

#set(EXECUTABLE_OUTPUT_PATH bin/)
set(COMPILE_FLAGS ${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})

add_executable(adapt ${SRC})

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
list(APPEND EXTERNAL_LIB ${MPI_LIBRARIES})

find_package(LAPACK)

find_path(HDF5_INCLUDE_DIR hdf5.h PATHS "${DEFAULT_HDF5_ROOT}/include")
if(NOT EXISTS "${HDF5_INCLUDE_DIR}")
  message(FATAL_ERROR "hdf5 include dir not found")
endif()

find_path(PARMETIS_INCLUDE_DIR parmetis.h PATHS "${DEFAULT_PARMETIS_ROOT}/include")
if(NOT EXISTS "${PARMETIS_INCLUDE_DIR}")
  message(FATAL_ERROR "parmetis include dir not found")
endif()

find_path(METIS_INCLUDE_DIR metis.h PATHS "${DEFAULT_METIS_ROOT}/include")
if(NOT EXISTS "${METIS_INCLUDE_DIR}")
  message(FATAL_ERROR "metis include dir not found")
endif()


find_library(HDF5_LIBRARY hdf5 PATHS "${DEFAULT_HDF5_ROOT}/lib")
if(NOT EXISTS "${HDF5_LIBRARY}")
  message(FATAL_ERROR "hdf5 library not found")
endif()

find_library(PARMETIS_LIBRARY parmetis PATHS "${DEFAULT_PARMETIS_ROOT}/lib")
if(NOT EXISTS "${PARMETIS_LIBRARY}")
  message(FATAL_ERROR "parmetis library not found")
endif()

find_library(METIS_LIBRARY metis PATHS "${DEFAULT_METIS_ROOT}/lib")
if(NOT EXISTS "${METIS_LIBRARY}")
  message(FATAL_ERROR "metis library not found")
endif()

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIR})

#list(APPEND EXTERNAL_LIB ${PARMETIS_LIBRARIES})
#list(APPEND EXTERNAL_LIB ${HDF5_LIBRARY})

target_link_libraries(adapt ${HDF5_LIBRARY})
target_link_libraries(adapt ${PARMETIS_LIBRARIES})
target_link_libraries(adapt ${MPI_LIBRARIES})
target_link_libraries(adapt ${LAPACK_LIBRARIES})
#====== Adding the external libraries here ======


