# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)

#set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_C_COMPILER icpc)

# set the project name and version
project(MIMIC VERSION 1.0)
# specify the C++ standard
include(./NAS_Dirk_work.cmake)
set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")
file(GLOB SRC CONFIGURE_DEPENDS "../../src/*.cpp" "main.cpp")
find_package(MPI REQUIRED)
#find_package(Boost COMPONENTS system filesystem REQUIRED)
#set(CMAKE_CXX_FLAGS "-Wl,-no_compact_unwind")

add_executable(partTetra ${SRC})

find_package(LAPACK)
#find_package(HDF5)
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

#if (VTK_FOUND)
#message("VTK is Found")
#endif()

#find_package(Boost COMPONENTS system filesystem REQUIRED)
set(Boost_INCLUDE_DIR /nobackupp19/dekelsch/boost_1_71_0/boost-install/include)
set(Boost_LIBRARY_DIR /nobackupp19/dekelsch/boost_1_71_0/boost-install/stage/lib)
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

#find_path(MPI_INCLUDE_DIR mpi.h PATHS "${DEFAULT_MPI_ROOT}/include")
#if(NOT EXISTS "${MPI_INCLUDE_DIR}")
#  message(FATAL_ERROR "mpi include dir not found")
#endif()

#find_path(CGNS_INCLUDE_DIR cgnslib.h PATHS "${DEFAULT_CGNS_ROOT}/include")
#if(NOT EXISTS "${CGNS_INCLUDE_DIR}")
#  message(FATAL_ERROR "cgns include dir not found")
#endif()





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

#find_library(MPI_LIBRARY mpi PATHS "${DEFAULT_MPI_ROOT}/lib")
#if(NOT EXISTS "${MPI_LIBRARY}")
#  message(FATAL_ERROR "mpi library not found")
#endif()

#find_library(CGNS_LIBRARY cgns PATHS "${DEFAULT_CGNS_ROOT}/lib")
#if(NOT EXISTS "${CGNS_LIBRARY}")
#  message(FATAL_ERROR "cgns library not found")
#endif()


set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${MMGS_INCLUDE_DIR})
include_directories(${MMG3D_INCLUDE_DIR})
include_directories(${PARMMG_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
#include_directories(${CGNS_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#list(APPEND EXTERNAL_LIB ${PARMETIS_LIBRARIES})
#list(APPEND EXTERNAL_LIB ${HDF5_LIBRARY})

#target_link_directories(partTetra INTERFACE /lib64)

target_link_libraries(partTetra ${VTK_LIBRARIES})
#target_link_libraries(partTetra ${VTK_LIBRARIES} "-shared /lib64/libm.so")
target_link_libraries(partTetra ${MMG_LIBRARY})
target_link_libraries(partTetra ${PARMMG_LIBRARY})
target_link_libraries(partTetra ${HDF5_LIBRARY})
target_link_libraries(partTetra ${PARMETIS_LIBRARIES})
target_link_libraries(partTetra ${MPI_CXX_LIBRARIES})
target_link_libraries(partTetra ${CGNS_LIBRARY})
target_link_libraries(partTetra ${LAPACK_LIBRARIES})
#target_link_libraries(partTetra ${Boost_LIBRARIES})
#target_link_libraries(partTetra ${Boost_LIBRARY_DIR})
#====== Adding the external libraries here ======



