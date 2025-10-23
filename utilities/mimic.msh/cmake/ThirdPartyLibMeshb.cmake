#FIND_LIBRARY(MMG_LIBRARY NAMES mmg)
#FIND_PATH(MMG_INCLUDE_DIR mmg.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(LIBMESHB_LIBRARY AND LIBMESHB_INCLUDE_DIR)
        SET(BUILD_LIBMESHB OFF)
ELSE()
        SET(BUILD_LIBMESHB ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_LIBMESHB
        "Build LibMeshb library from ThirdParty." ${BUILD_LIBMESHB})

IF (THIRDPARTY_BUILD_LIBMESHB)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build libMeshb.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                libMeshb
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
		GIT_REPOSITORY https://github.com/LoicMarechal/libMeshb.git
    		URL_MD5 dbabd502ff1c7946ed903f2cbf9355a0
                STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/libMeshb
            	BINARY_DIR ${TPBUILD}/libMeshb
            	TMP_DIR ${TPBUILD}/libMeshb-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/libMeshb
                )
        THIRDPARTY_LIBRARY(LIBMESHB_LIBRARY STATIC Meshb.7
            DESCRIPTION "LIBMESHB library")

        SET(LIBMESHB_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "LibMeshb include directory" FORCE)
        MESSAGE(STATUS "Build LIBMESHB: ${LIBMESHB_LIBRARY}")
        SET(MMG_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
#ADD_DEPENDENCIES(thirdparty MMG)
MARK_AS_ADVANCED(LIBMESHB_LIBRARY)
MARK_AS_ADVANCED(LIBMESHB_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${LIBMESHB_INCLUDE_DIR})


