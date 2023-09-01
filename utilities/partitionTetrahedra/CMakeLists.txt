# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)

# set the project name and version
project(MIMIC VERSION 1.0)


# specify the C++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

include(../../platforms_cmake/MacOSX_Dirk_work_serialHDF5.cmake)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wno-dev)

#set(LIBRARY_OUTPUT_PATH lib/)

set(SRC ../../src/NekFace.cpp ../../src/tinyxmlparser.cpp ../../src/tinyxmlerror.cpp ../../src/tinystr.cpp ../../src/tinyxml.cpp ../../src/adapt_parops.cpp ../../src/adapt_bltopology.cpp ../../src/adapt_boundary.cpp ../../src/adapt_output.cpp ../../src/adapt_compute.cpp ../../src/adapt_schedule.cpp ../../src/adapt_operations.cpp ../../src/hex2tet.cpp ../../src/adapt_geometry.cpp ../../src/adapt_math.cpp ../../src/adapt_io.cpp ../../src/adapt_recongrad.cpp ../../src/adapt_topology.cpp ../../src/adapt_partition.cpp ../../src/adapt_prismaticlayer.cpp ../../src/adapt_redistribute.cpp ../../src/adapt_prismaticlayer.cpp ../../src/adapt_parstate.cpp ../../src/adapt_distri_parstate.cpp main.cpp)

#set(EXECUTABLE_OUTPUT_PATH bin/)
set(COMPILE_FLAGS ${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_FLAGS "-Wl,-no_compact_unwind")
add_executable(us3d2cgns ${SRC})

#find_package(MPI REQUIRED)
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
#list(APPEND EXTERNAL_LIB ${MPI_LIBRARIES})

find_package(LAPACK)

#find_path(XML_INCLUDE_DIR PATHS "${DEFAULT_XML_ROOT}")
#if(NOT EXISTS "${XML_INCLUDE_DIR}")
#  message(FATAL_ERROR "tinyxml include dir not found")
#endif()

set(Boost_INCLUDE_DIR /Users/dekelsch/Software/boost_1_71_0)
set(Boost_LIBRARY_DIR /Users/dekelsch/Software/boost_1_71_0/stage/lib)
find_package(Boost COMPONENTS system filesystem REQUIRED)

find_path(MMGS_INCLUDE_DIR /mmg/mmgs/libmmgs.h PATHS "${DEFAULT_MMG_ROOT}/include")
if(NOT EXISTS "${MMGS_INCLUDE_DIR}")
  message(FATAL_ERROR "libmmgs include dir not found")
endif()

find_path(MMG3D_INCLUDE_DIR /mmg/mmg3d/libmmg3d.h PATHS "${DEFAULT_MMG_ROOT}/include")
if(NOT EXISTS "${MMG3D_INCLUDE_DIR}")
  message(FATAL_ERROR "libmmg3d include dir not found")
endif()

find_path(PARMMG_INCLUDE_DIR /parmmg/libparmmg.h PATHS "${DEFAULT_PARMMG_ROOT}/include")
if(NOT EXISTS "${PARMMG_INCLUDE_DIR}")
  message(FATAL_ERROR "libparmmg include dir not found")
endif()

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

find_path(MPI_INCLUDE_DIR mpi.h PATHS "${DEFAULT_MPI_ROOT}/include")
if(NOT EXISTS "${MPI_INCLUDE_DIR}")
  message(FATAL_ERROR "mpi include dir not found")
endif()

find_path(CGNS_INCLUDE_DIR cgnslib.h PATHS "${DEFAULT_CGNS_ROOT}/include")
if(NOT EXISTS "${CGNS_INCLUDE_DIR}")
  message(FATAL_ERROR "cgns include dir not found")
endif()

find_library(MMG_LIBRARY mmg PATHS "${DEFAULT_MMG_ROOT}/lib")
if(NOT EXISTS "${MMG_LIBRARY}")
  message(FATAL_ERROR "mmg library not found")
endif()

find_library(PARMMG_LIBRARY parmmg PATHS "${DEFAULT_PARMMG_ROOT}/lib")
if(NOT EXISTS "${PARMMG_LIBRARY}")
  message(FATAL_ERROR "parmmg library not found")
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

find_library(MPI_LIBRARY mpi PATHS "${DEFAULT_MPI_ROOT}/lib")
if(NOT EXISTS "${MPI_LIBRARY}")
  message(FATAL_ERROR "mpi library not found")
endif()

find_library(CGNS_LIBRARY cgns PATHS "${DEFAULT_CGNS_ROOT}/lib")
if(NOT EXISTS "${CGNS_LIBRARY}")
  message(FATAL_ERROR "cgns library not found")
endif()


set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${MMGS_INCLUDE_DIR})
include_directories(${MMG3D_INCLUDE_DIR})
include_directories(${PARMMG_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_DIR})
include_directories(${CGNS_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#list(APPEND EXTERNAL_LIB ${PARMETIS_LIBRARIES})
#list(APPEND EXTERNAL_LIB ${HDF5_LIBRARY})

target_link_libraries(us3d2cgns ${MMG_LIBRARY})
target_link_libraries(us3d2cgns ${PARMMG_LIBRARY})
target_link_libraries(us3d2cgns ${HDF5_LIBRARY})
target_link_libraries(us3d2cgns ${PARMETIS_LIBRARIES})
target_link_libraries(us3d2cgns ${MPI_LIBRARY})
target_link_libraries(us3d2cgns ${CGNS_LIBRARY})
target_link_libraries(us3d2cgns ${LAPACK_LIBRARIES})
target_link_libraries(us3d2cgns ${Boost_LIBRARIES} )

#====== Adding the external libraries here ======


