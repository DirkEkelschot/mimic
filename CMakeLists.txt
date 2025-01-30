# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)
set(CMAKE_CXX_STANDARD 17)
# set the project name and version
project(mimic VERSION 1.0)
# specify the C++ standard
include(platforms_cmake/MacOSX_Dirk_work.cmake)
include(ExternalProject)
INCLUDE(CMakeDependentOption)


find_package(LAPACK)
MESSAGE("${MPI_CXX_LIBRARIES} annnad ${MPI_CXX_INCLUDE_PATH}")
find_package(VTK 
	COMPONENTS 
	CommonCore
        CommonDataModel
        IOLegacy
        IOXML
        ParallelCore
        IOParallelXML
	REQUIRED
)

if (VTK_FOUND)
message("VTK is Found")
endif()

find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(MPI REQUIRED)

INCLUDE (cmake/adaptCommon.cmake)
#INCLUDE (cmake/ThirdPartyTetGen.cmake)
INCLUDE (cmake/ThirdPartyGSL.cmake)
#INCLUDE (cmake/ThirdPartyLibMeshb.cmake)
INCLUDE (cmake/ThirdPartyZLIB.cmake)
INCLUDE (cmake/ThirdPartyHDF5.cmake)
INCLUDE (cmake/ThirdPartyCGAL.cmake)
INCLUDE (cmake/ThirdPartyMMG.cmake)
INCLUDE (cmake/ThirdPartyParMMG.cmake)
#find_path(TETGEN_INCLUDE_DIR tetgen.h PATHS "${DEFAULT_TETGEN_ROOT}")
#if(NOT EXISTS "${TETGEN_INCLUDE_DIR}")
#  message(FATAL_ERROR "tetgen include dir not found")
#endif()

#message(${TETGEN_LIBRARY})
message("HAHAHA",${DEFAULT_VORO_ROOT})

file(GLOB SRC CONFIGURE_DEPENDS "src/*.cpp" "main.cpp")
set(CMAKE_CXX_FLAGS "-Wl,-no_compact_unwind")
add_executable(mimic ${SRC})



find_path(PARMETIS_INCLUDE_DIR parmetis.h PATHS "${DEFAULT_PARMETIS_ROOT}/include")
if(NOT EXISTS "${PARMETIS_INCLUDE_DIR}")
  message(FATAL_ERROR "parmetis include dir not found")
endif()

find_path(METIS_INCLUDE_DIR metis.h PATHS "${DEFAULT_METIS_ROOT}/include")
if(NOT EXISTS "${METIS_INCLUDE_DIR}")
  message(FATAL_ERROR "metis include dir not found")
endif()


#find_path(VORO_INCLUDE_DIR voro++.hh PATHS "${DEFAULT_VORO_ROOT}/include/voro++")
#if(NOT EXISTS "${VORO_INCLUDE_DIR}")
#  message(FATAL_ERROR "voro++ include dir not found")
#endif()

find_library(PARMETIS_LIBRARY parmetis PATHS "${DEFAULT_PARMETIS_ROOT}/lib")
if(NOT EXISTS "${PARMETIS_LIBRARY}")
  message(FATAL_ERROR "parmetis library not found")
endif()

find_library(METIS_LIBRARY metis PATHS "${DEFAULT_METIS_ROOT}/lib")
if(NOT EXISTS "${METIS_LIBRARY}")
  message(FATAL_ERROR "metis library not found")
endif()


#find_path(HDF5_INCLUDE_DIR hdf5.h PATHS "${DEFAULT_HDF5_ROOT}/include")
#if(NOT EXISTS "${HDF5_INCLUDE_DIR}")
#  message(FATAL_ERROR "hdf5 include dir not found")
#endif()




find_package(CGAL REQUIRED)


#find_library(HDF5_LIBRARY hdf5 PATHS "${DEFAULT_HDF5_ROOT}/lib")
#if(NOT EXISTS "${HDF5_LIBRARY}")
#  message(FATAL_ERROR "hdf5 library not found")
#endif()

#message(${DEFAULT_VORO_ROOT}/lib)
#find_path(VORO_LIBRARY libvoro++.a PATHS "${DEFAULT_VORO_ROOT}/lib")
#if(NOT EXISTS "${VORO_LIBRARY}")
#  message(FATAL_ERROR "voro++ lib dir not found#")
#endif()



include_directories(${MPI_CXX_INCLUDE_DIR})

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
#include_directories(${MMGS_INCLUDE_DIR})
include_directories(${MMG3D_INCLUDE_DIR})
#include_directories(${TETGEN_INCLUDE_DIRS})
include_directories(${PARMMG_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
#include_directories(${CGNS_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
#include_directories(${VORO_INCLUDE_DIR})
#include_directories(${CGAL_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
include_directories(${GSL_INCLUDE_DIRS})

#list(APPEND EXTERNAL_LIB ${PARMETIS_LIBRARIES})
#list(APPEND EXTERNAL_LIB ${HDF5_LIBRARY})
target_link_libraries(mimic ${VTK_LIBRARIES})
#target_link_libraries(mimic ${TETGEN_LIBRARY})
target_link_libraries(mimic ${MMG_LIBRARY})
target_link_libraries(mimic ${GSL_LIBRARY})
#target_link_libraries(partTetra ${CMAKE_SOURCE_DIR}/build/ThirdParty/ParMmg/Mmg-prefix/src/Mmg-build/lib)
#target_link_libraries(mimic ${CGAL_LIBRARY})
target_link_libraries(mimic ${LIBMESHB_LIBRARY})
target_link_libraries(mimic ${PARMMG_LIBRARY})
target_link_libraries(mimic ${HDF5_LIBRARIES})
target_link_libraries(mimic ${ZLIB_LIBRARY})
#target_link_libraries(partTetra ${HDF5_LIBRARIES_SHARED})
#target_link_libraries(partTetra ${HDF5_LIBRARIES_STATIC})
target_link_libraries(mimic ${PARMETIS_LIBRARIES})
target_link_libraries(mimic ${MPI_CXX_LIBRARIES})
#target_link_libraries(readLibMeshb ${CGNS_LIBRARY})
target_link_libraries(mimic ${LAPACK_LIBRARIES})
target_link_libraries(mimic ${Boost_LIBRARIES} )
#target_link_libraries(mimic ${VORO_LIBRARY} )
target_link_libraries(mimic CGAL::CGAL)
#====== Adding the external libraries here ======
