

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(GSL_LIBRARY AND GSL_INCLUDE_DIR)
        SET(BUILD_GSL OFF)
ELSE()
        SET(BUILD_GSL ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_GSL
        "Build TetGen library from ThirdParty." ${BUILD_GSL})

IF (THIRDPARTY_BUILD_GSL)
        INCLUDE(ExternalProject)
#        UNSET(PATCH CACHE)
#        FIND_PROGRAM(PATCH patch)
#        IF(NOT PATCH)
#                MESSAGE(FATAL_ERROR
#                "'patch' tool for modifying files not found. Cannot build gsl.")
#        ENDIF()
#        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                gsl-2.8
		URL https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
		URL_MD5 182ec03204f164e67238c9116591a37d
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/gsl-2.8
                BINARY_DIR ${TPBUILD}/gsl-2.8
                TMP_DIR ${TPBUILD}/gsl-2.8-tmp
		CONFIGURE_COMMAND ${TPSRC}/gsl-2.8/./configure --prefix=${TPDIST} CC=mpicc CXX=mpic++
		BUILD_COMMAND make all
#		INSTALL_COMMAND make install
                ${TPSRC}/gsl-2.8
                )
        THIRDPARTY_LIBRARY(GSL_LIBRARY STATIC gsl
            DESCRIPTION "GSL library")
	
	SET(GSL_DIR ${TPDIST} CACHE FILEPATH
            "GSL DIR directory" FORCE)
        SET(GSL_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "GSL include directory" FORCE)
	SET(GSL_INCLUDE_LIB ${TPDIST}/lib CACHE FILEPATH
            "GSL lib directory" FORCE)
        MESSAGE(STATUS "Build GSL: ${GSL_LIBRARY}")
        SET(GSL_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
MARK_AS_ADVANCED(GSL_LIBRARY)
MARK_AS_ADVANCED(GSL_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})


