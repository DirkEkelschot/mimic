# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)

#set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_STANDARD 17)



# set the project name and version
project(MIMIC VERSION 1.0)
# specify the C++ standard
include(./machine.cmake)
set(CMAKE_C_COMPILER ${DEFAULT_MPI_ROOT}/bin/mpicc)
set(CMAKE_CXX_COMPILER ${DEFAULT_MPI_ROOT}/bin/mpicxx)
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")


file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")


	



if(APPLE)
    message("Building on macOS or iOS")
    set(LIB_EXT ".dylib" CACHE STRING "Dynamic library extension")
    # assuming running on Apple Silicon: arm64
    set(CMAKE_OSX_ARCHITECTURES arm64)
    # if building on Intel, set the architecture to x86_64
    # set(CMAKE_OSX_ARCHITECTURES arm64)
ELSEIF(UNIX)
    message("Building on Linux")
    set(LIB_EXT ".so" CACHE STRING "Dynamic library extension")
    set(CMAKE_OSX_ARCHITECTURES x86_64)
ENDIF()


set(_user_dep_vars
  DEFAULT_HDF5_ROOT
  DEFAULT_MMG_ROOT
  DEFAULT_PARMMG_ROOT
  DEFAULT_PARMETIS_ROOT
  VTK_DIR
  DEFAULT_BOOST_ROOT
)

foreach(var IN LISTS _user_dep_vars)
  if(DEFINED ${var})
    message(STATUS "${var} is set to '${${var}}'")
  else()
    message(STATUS "${var} is NOT set")
  endif()
endforeach()



INCLUDE (cmake/adaptCommon.cmake)




find_package(LAPACK)

set(Boost_INCLUDE_DIR ${DEFAULT_BOOST_ROOT}/include)
set(Boost_LIBRARY_DIR ${DEFAULT_BOOST_ROOT}/lib)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#=================================================================================================



set(MPI_C ${DEFAULT_MPI_ROOT}/bin/mpicc)
set(MPI_CXX ${DEFAULT_MPI_ROOT}/bin/mpicxx)

find_path(MPI_INCLUDE_DIR mpi.h PATHS "${DEFAULT_MPI_ROOT}/include")
if(NOT EXISTS "${MPI_INCLUDE_DIR}")
  message(FATAL_ERROR "mpi include dir not found")
endif()

find_library(MPI_LIBRARY mpi PATHS "${DEFAULT_MPI_ROOT}/lib")
if(NOT EXISTS "${MPI_LIBRARY}")
  message(FATAL_ERROR "mpi library not found")
endif()
include_directories(${MPI_INCLUDE_DIR})




#=================================================================================================


IF(DEFAULT_PARMETIS_ROOT AND DEFAULT_METIS_ROOT)

  find_path(PARMETIS_INCLUDE_DIR parmetis.h PATHS "${DEFAULT_PARMETIS_ROOT}/include")
  if(NOT EXISTS "${PARMETIS_INCLUDE_DIR}")
    message(FATAL_ERROR "parmetis include dir not found")
  endif()

  find_path(METIS_INCLUDE_DIR metis.h PATHS "${DEFAULT_METIS_ROOT}/include")
  if(NOT EXISTS "${METIS_INCLUDE_DIR}")
    message(FATAL_ERROR "metis include dir not found")
  endif()


  find_library(PARMETIS_LIBRARY parmetis PATHS "${DEFAULT_PARMETIS_ROOT}/lib")
  if(NOT EXISTS "${PARMETIS_LIBRARY}")
    message(FATAL_ERROR "parmetis library not found")
  endif()

  find_library(METIS_LIBRARY metis PATHS "${DEFAULT_METIS_ROOT}/lib")
  if(NOT EXISTS "${METIS_LIBRARY}")
    message(FATAL_ERROR "metis library not found")
  endif()

ELSEIF(NOT DEFAULT_PARMETIS_ROOT OR NOT DEFAULT_METIS_ROOT)
  MESSAGE(FATAL_ERROR [=[
Configuration failed!
Install METIS and PARMETIS and set DEFAULT_METIS_ROOT and DEFAULT_PARMETIS_ROOT equal to the installing paths.
For more information on how to install METIS and PARMETIS, see the README.md file.
Modify the paths listed below in machine.cmake:
 - DEFAULT_METIS_ROOT=/path/to/metis/
 - DEFAULT_PARMETIS_ROOT=/path/to/parmetis/
Make sure you have appropriate license arrangements in case you are using this software for commercial purposes.
Please install these and try again.
]=])
ENDIF()


set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})

include_directories(${PARMETIS_INCLUDE_DIRS})


#=================================================================================================


IF(DEFAULT_HDF5_ROOT)

  if(HDF5_VERSION VERSION_LESS "1.10.0")
    message(WARNING "HDF5 version ${HDF5_VERSION} is too old. Please use HDF5 >= 1.10.0 for proper support.")
  endif()
  
  find_path(HDF5_INCLUDE_DIR hdf5.h PATHS "${DEFAULT_HDF5_ROOT}/include")
  if(NOT EXISTS "${HDF5_INCLUDE_DIR}")
    message(FATAL_ERROR "hdf5 include dir not found")
  endif()

  find_library(HDF5_LIBRARY hdf5 PATHS "${DEFAULT_HDF5_ROOT}/lib")
  if(NOT EXISTS "${HDF5_LIBRARY}")
    message(FATAL_ERROR "hdf5 library not found")
  endif()

  find_library(HDF5_HL_LIBRARY hdf5_hl PATHS "${DEFAULT_HDF5_ROOT}/lib")
  if(NOT EXISTS "${HDF5_LIBRARY}")
    message(FATAL_ERROR "hdf5_hl library not found")
  endif()

  SET(HDF5_DIR ${DEFAULT_HDF5_ROOT})
  SET(HDF5_INCLUDE_DIRS ${HDF5_LIBRARY})
  SET(HDF5_LIBRARIES ${HDF5_INCLUDE_DIR})
  include_directories(${HDF5_INCLUDE_DIR})

ENDIF()




#=================================================================================================


IF(DEFAULT_MMG_ROOT)
  
  find_path(MMG3D_INCLUDE_DIR /mmg/mmg3d/libmmg3d_private.h PATHS "${DEFAULT_MMG_ROOT}/include")
  if(NOT EXISTS "${MMG3D_INCLUDE_DIR}")
    message(FATAL_ERROR "libmmg3d include dir not found")
  endif()

  find_path(MMGS_INCLUDE_DIR /mmg/mmgs/libmmgs.h PATHS "${DEFAULT_MMG_ROOT}/include")
  if(NOT EXISTS "${MMGS_INCLUDE_DIR}")
    message(FATAL_ERROR "libmmgs include dir not found")
  endif()

  find_library(MMG_LIBRARY mmg3d PATHS "${DEFAULT_MMG_ROOT}/lib")
  if(NOT EXISTS "${MMG_LIBRARY}")
    message(FATAL_ERROR "mmg3d library not found")
  endif()

  include_directories(${MMGS_INCLUDE_DIR})
  include_directories(${MMG3D_INCLUDE_DIR})

ELSE()



ENDIF()




#=================================================================================================


IF(DEFAULT_PARMMG_ROOT)
  
  find_path(PARMMG_INCLUDE_DIR parmmg.h PATHS "${DEFAULT_PARMMG_ROOT}/include")
  if(NOT EXISTS "${PARMMG_INCLUDE_DIR}")
    message(FATAL_ERROR "parmmg include dir not found")
  endif()

  find_library(PARMMG_LIBRARY parmmg PATHS "${DEFAULT_PARMMG_ROOT}/lib")
  if(NOT EXISTS "${PARMMG_LIBRARY}")
    message(FATAL_ERROR "parmmg library not found")
  endif()

ELSE()

  IF(NOT DEFAULT_HDF5_ROOT)
    MESSAGE("Installing HDF5 in the ThirdPary directory...")
    MESSAGE("HDF5 requires ZLIB as well and ZLIB is installed in the ThirdParty directory.")

    MESSAGE("Loading cmake/ThirdPartyZLIB.cmake")
    INCLUDE (cmake/ThirdPartyZLIB.cmake)

    MESSAGE("Loading cmake/ThirdPartyHDF5.cmake")
    INCLUDE (cmake/ThirdPartyHDF5.cmake)

    include_directories(${HDF5_INCLUDE_DIR})

    IF(NOT DEFAULT_MMG_ROOT)
      MESSAGE("Loading cmake/ThirdPartyMMG.cmake")
      INCLUDE (cmake/ThirdPartyMMG.cmake)

      include_directories(${MMGS_INCLUDE_DIR})
      include_directories(${MMG3D_INCLUDE_DIR})

    ENDIF()
  
    MESSAGE("Loading cmake/ThirdPartyParMMG.cmake")
    INCLUDE (cmake/ThirdPartyParMMG.cmake)
  ENDIF()

  IF(DEFAULT_HDF5_ROOT)
  
    MESSAGE("HDF5 found in user specified path:" ${HDF5_DIR})
    IF(NOT DEFAULT_MMG_ROOT)
      MESSAGE("Loading cmake/ThirdPartyMMG.cmake")
      INCLUDE (cmake/ThirdPartyMMG.cmake)

      include_directories(${MMGS_INCLUDE_DIR})
      include_directories(${MMG3D_INCLUDE_DIR})

    ENDIF()
  
    MESSAGE("Loading cmake/ThirdPartyParMMG.cmake")
    INCLUDE (cmake/ThirdPartyParMMG.cmake)
    
  
  ENDIF()



ENDIF()

include_directories(${PARMMG_INCLUDE_DIR})

#=================================================================================================



IF(VTK_DIR)


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


  MESSAGE("VTK found in user specified path:" ${VTK_LIBRARY})

ELSE()

  MESSAGE("Loading cmake/ThirdPartyVTK.cmake")
  INCLUDE (cmake/ThirdPartyVTK.cmake)

	
set(VTK_LIBRARIES
  ${VTK_LIBRARY}/libvtksys-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkCommonCore-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkCommonDataModel-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOLegacy-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOXML-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkParallelCore-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOParallelXML-9.4${LIB_EXT}
)

include_directories(${VTK_INCLUDE_DIRS})

ENDIF()

OPTION(BUILD_MIMIC_US3D
        "Build MIMIC for US3D." OFF)
OPTION(BUILD_MIMIC_HYPERSOLVE
        "Build MIMIC for HyperSolve." OFF)

#add_subdirectory(unit_tests)

IF(BUILD_MIMIC_US3D)
	add_subdirectory(mimic.us3d)
ENDIF()
IF(BUILD_MIMIC_HYPERSOLVE)
	add_subdirectory(mimic.hs)
ENDIF()



OPTION(BUILD_UTILITIES
        "Build MIMIC utilities." OFF)

IF(BUILD_UTILITIES)
	add_subdirectory(utilities)
ENDIF()




#=================================================================================================


















