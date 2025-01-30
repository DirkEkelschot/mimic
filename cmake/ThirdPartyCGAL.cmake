FIND_LIBRARY(CGAL_LIBRARY NAMES cgal)
FIND_PATH(CGAL_INCLUDE_DIR cgal.h)
IF(CGAL_LIBRARY AND CGAL_INCLUDE_DIR)
        SET(BUILD_CGAL OFF)
ELSE()
        SET(BUILD_CGAL ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_CGAL
        "Build CGAL library from ThirdParty." ${BUILD_CGAL})

IF (THIRDPARTY_BUILD_CGAL)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build CGAL.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                CGAL-5.6
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL https://github.com/CGAL/cgal/releases/download/v5.6/CGAL-5.6.zip
                URL_MD5 6d1d067b88e20f7080d07d5108b4c772
                STAMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/stamp
                DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/Thirdparty
                SOURCE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/CGAL-5.6
                BINARY_DIR ${CMAKE_BINARY_DIR}/Thirdparty/CGAL-5.6
                TMP_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/CGAL-5.6-tmp
                INSTALL_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist
                CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Thirdparty/dist
                -DGMP_INCLUDE_DIR=/Users/dekelsch/brew-4.1.24/Cellar/gmp/6.3.0/include
		-DGMP_LIBRARIES=/Users/dekelsch/brew-4.1.24/Cellar/gmp/6.3.0/lib
		${CMAKE_BINARY_DIR}/Thirdparty/CGAL-5.6
                )
        THIRDPARTY_LIBRARY(CGAL_LIBRARY STATIC cgal
            DESCRIPTION "CGAL library")
        SET(CGAL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/include CACHE FILEPATH
            "CGAL include" FORCE)
        ADD_DEFINITIONS(-DCGAL_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build CGAL: ${CGAL_LIBRARY}")
        SET(CGAL_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
message(${CGAL_INCLUDE_DIR})
MARK_AS_ADVANCED(CGAL_LIBRARY)
MARK_AS_ADVANCED(CGAL_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
