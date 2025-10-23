#FIND_LIBRARY(HDF5_LIBRARY NAMES hdf5)
#FIND_PATH(TETGEN_INCLUDE_DIR hdf5.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(HDF5_LIBRARY AND HDF5_INCLUDE_DIR)
        SET(BUILD_HDF5 OFF)
ELSE()
        SET(BUILD_HDF5 ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_HDF5
        "Build TetGen library from ThirdParty." ${BUILD_HDF5})

IF (THIRDPARTY_BUILD_HDF5)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                hdf5-1_14_3
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                #URL https://www.nektar.info/thirdparty/hdf5-1.12.3.tar.bz2
                #URL_MD5 5d609bf2a74f980aa42dbe61de452185
                URL https://github.com/HDFGroup/hdf5/releases/download/hdf5-1_14_3/hdf5-1_14_3.tar.gz
		URL_MD5 d862cc7534526549ecc544df89286f0b
		STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/hdf5-1_14_3
            	BINARY_DIR ${TPBUILD}/hdf5-1_14_3 
            	TMP_DIR ${TPBUILD}/hdf5-1_14_3-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DHDF5_ENABLE_PARALLEL=ON
		-DZLIB_DIR=${TPDIST}
                ${TPSRC}/hdf5-1_14_3
                )
        THIRDPARTY_LIBRARY(HDF5_LIBRARY STATIC hdf5
            DESCRIPTION "HDF5 library shared")
	#THIRDPARTY_LIBRARY(HDF5_LIBRARIES_STATIC STATIC hdf5
        #    DESCRIPTION "HDF5 library static")
        SET(HDF5_INCLUDE_DIRS ${TPDIST}/include CACHE FILEPATH
            "HDF5 include directory" FORCE)
        MESSAGE(STATUS "Build HDF5 shared: ${HDF5_LIBRARY}")
        #MESSAGE(STATUS "Build HDF5 static: ${HDF5_LIBRARIES_STATIC}")
	SET(HDF5_CONFIG_INCLUDE_DIR ${TPINC})

ENDIF()
#SET(HDF5_LIB ${TPDIST}/lib)
SET(HDF5_LIBRARIES ${TPDIST}/lib/libhdf5.a)
MARK_AS_ADVANCED(HDF5_LIBRARIES)
#MARK_AS_ADVANCED(HDF5_LIBRARIES_STATIC)
MARK_AS_ADVANCED(HDF5_INCLUDE_DIRS)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})


