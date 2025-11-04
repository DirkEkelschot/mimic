

#FIND_LIBRARY(TETGEN_LIBRARY NAMES tet)
#FIND_PATH(TETGEN_INCLUDE_DIR tetgen.h)


SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(TETGEN_LIBRARY AND TETGEN_INCLUDE_DIR)
        SET(BUILD_TETGEN OFF)
ELSE()
        SET(BUILD_TETGEN ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_TETGEN
        "Build TetGen library from ThirdParty." ${BUILD_TETGEN})

IF (THIRDPARTY_BUILD_TETGEN)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                tetgen-1.5
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL https://www.nektar.info/thirdparty/tetgen-1.5.zip
                URL_MD5 6d62e63f9b1e7a8ce53d5bc87e6a0a09
                STAMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/stamp
                DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/Thirdparty
                SOURCE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/tetgen-1.5
                BINARY_DIR ${CMAKE_BINARY_DIR}/Thirdparty/tetgen-1.5
                TMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/tetgen-1.5-tmp
                INSTALL_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist
                CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Thirdparty/dist
                ${CMAKE_BINARY_DIR}/Thirdparty/tetgen-1.5
                )
        THIRDPARTY_LIBRARY(TETGEN_LIBRARY STATIC tetgen
            DESCRIPTION "Tetgen library")
        SET(TETGEN_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/include CACHE FILEPATH
            "TetGen include" FORCE)
        ADD_DEFINITIONS(-DTETGEN_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build TetGen: ${TETGEN_LIBRARY}")
        SET(TETGEN_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
message(${TETGEN_INCLUDE_DIR})
MARK_AS_ADVANCED(TETGEN_LIBRARY)
MARK_AS_ADVANCED(TETGEN_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${TETGEN_INCLUDE_DIR})
