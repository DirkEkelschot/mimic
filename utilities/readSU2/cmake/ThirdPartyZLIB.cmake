SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

FIND_LIBRARY(ZLIB_LIBRARY NAMES zlib)
FIND_PATH(ZLIB_INCLUDE_DIR zlib.h)
IF(ZLIB_LIBRARY AND ZLIB_INCLUDE_DIR)
        SET(BUILD_ZLIB OFF)
ELSE()
        SET(BUILD_ZLIB ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_ZLIB
        "Build ZLIB library from ThirdParty." ${BUILD_ZLIB})

IF (THIRDPARTY_BUILD_ZLIB)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build ZLIB.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                zlib-1.3.1
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
                URL_MD5 9855b6d802d7fe5b7bd5b196a2271655
                STAMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/stamp
                DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/Thirdparty
                SOURCE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/zlib-1.3.1
                BINARY_DIR ${CMAKE_BINARY_DIR}/Thirdparty/zlib-1.3.1
                TMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/zlib-1.3.1-tmp
                INSTALL_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist
                CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Thirdparty/dist
                ${CMAKE_BINARY_DIR}/Thirdparty/zlib-1.3.1
                )
        #THIRDPARTY_LIBRARY(ZLIB_LIBRARY STATIC z
        #    DESCRIPTION "ZLIB library")
        SET(ZLIB_LIBRARY ${TPDIST}/lib/libz.a CACHE FILEPATH
            "ZLIB lib" FORCE)
	SET(ZLIB_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/include CACHE FILEPATH
            "ZLIB include" FORCE)
        ADD_DEFINITIONS(-DZLIB_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build ZLIB: ${ZLIB_LIBRARY}")
        SET(ZLIB_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
MESSAGE("----fiun--> ${ZLIB_LIBRARY}------<")
message(${ZLIB_INCLUDE_DIR})
MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
