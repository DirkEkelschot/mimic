

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(GMP_LIBRARY AND GMP_INCLUDE_DIR)
        SET(BUILD_GMP OFF)
ELSE()
        SET(BUILD_GMP ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_GMP
        "Build TetGen library from ThirdParty." ${BUILD_GMP})

IF (THIRDPARTY_BUILD_GMP)
        INCLUDE(ExternalProject)
#        UNSET(PATCH CACHE)
#        FIND_PROGRAM(PATCH patch)
#        IF(NOT PATCH)
#                MESSAGE(FATAL_ERROR
#                "'patch' tool for modifying files not found. Cannot build gmp.")
#        ENDIF()
#        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                gmp-6.3.0
		URL https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/gmp-6.3.0
                BINARY_DIR ${TPBUILD}/gmp-6.3.0
                TMP_DIR ${TPBUILD}/gmp-6.3.0-tmp
		CONFIGURE_COMMAND ${TPSRC}/gmp-6.3.0/./configure --prefix=${TPDIST} CC=mpicc CXX=mpic++
		BUILD_COMMAND make all
#		INSTALL_COMMAND make install
                ${TPSRC}/gmp-6.3.0
                )
        THIRDPARTY_LIBRARY(GMP_LIBRARY STATIC gmp
            DESCRIPTION "GMP library")
	
	SET(GMP_DIR ${TPDIST} CACHE FILEPATH
            "GMP DIR directory" FORCE)
        SET(GMP_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "GMP include directory" FORCE)
	SET(GMP_INCLUDE_LIB ${TPDIST}/lib CACHE FILEPATH
            "GMP lib directory" FORCE)
        MESSAGE(STATUS "Build GMP: ${GMP_LIBRARY}")
        SET(GMP_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
SET(GMP_LIBRARY ${TPDIST}/lib/libgmp.a)
MARK_AS_ADVANCED(GMP_LIBRARY)
MARK_AS_ADVANCED(GMP_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})


