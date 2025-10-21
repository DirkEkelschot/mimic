FIND_LIBRARY(PARMETIS_LIBRARY NAMES parmetis)
FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h)
IF(PARMETIS_LIBRARY AND PARMETIS_INCLUDE_DIR)
        SET(BUILD_PARMETIS OFF)
ELSE()
        SET(BUILD_PARMETIS ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_PARMETIS
        "Build Parmetis library from ThirdParty." ${BUILD_TETGEN})

IF (THIRDPARTY_BUILD_PARMETIS)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Parmetis.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                parmetis-4.0.3
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
                URL_MD5 f69c479586bf6bb7aff6a9bc0c739628
                STAMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/stamp
                DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/Thirdparty
                SOURCE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/parmetis-4.0.3
                BINARY_DIR ${CMAKE_BINARY_DIR}/Thirdparty/parmetis-4.0.3
                TMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/parmetis-4.0.3-tmp
                INSTALL_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist
                CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Thirdparty/dist
                ${CMAKE_BINARY_DIR}/Thirdparty/parmetis-4.0.3
                )
        THIRDPARTY_LIBRARY(PARMETIS_LIBRARY STATIC parmetis
            DESCRIPTION "Parmetis library")
        SET(PARMETIS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/include CACHE FILEPATH
            "Parmetis include" FORCE)
        ADD_DEFINITIONS(-DTETGEN_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build Parmetis: ${PARMETIS_LIBRARY}")
        SET(PARMETIS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
message(${PARMETIS_INCLUDE_DIR})
MARK_AS_ADVANCED(PARMETIS_LIBRARY)
MARK_AS_ADVANCED(PARMETIS_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIR})
