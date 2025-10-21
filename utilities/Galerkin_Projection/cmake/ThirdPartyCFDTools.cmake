SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
message("CMAKE_SOURCE_DIR=" ${CMAKE_SOURCE_DIR}) 
FIND_LIBRARY(CFDTOOLS_LIBRARY NAMES CFDTOOLS)
FIND_PATH(CFDTOOLS_INCLUDE_DIR CFDTOOLS.h)
IF(CFDTOOLS_LIBRARY AND CFDTOOLS_INCLUDE_DIR)
        SET(BUILD_CFDTOOLS OFF)
ELSE()
        SET(BUILD_CFDTOOLS ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_CFDTOOLS
        "Build CFDTOOLS library from ThirdParty." ${BUILD_CFDTOOLS})

IF (THIRDPARTY_BUILD_CFDTOOLS)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build CFDTOOLS.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                cfdtools
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                GIT_REPOSITORY https://www.github.com/nasa/cfdtools.git
                GIT_TAG develop
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/cfdtools
                BINARY_DIR ${TPBUILD}/cfdtools
                TMP_DIR ${TPBUILD}/cfdtools-tmp
                INSTALL_DIR ${TPDIST}
                CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                ${TPSRC}/cfdtools
                )
        SET(CFDTOOLS_LIBRARY ${TPDIST}/lib/libcfdtools_kdtree.a CACHE FILEPATH
            "CFDTOOLS lib" FORCE)
	SET(CFDTOOLS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/Thirdparty/dist/include/ CACHE FILEPATH
            "CFDTOOLS include" FORCE)
        ADD_DEFINITIONS(-DCFDTOOLS_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build CFDTOOLS: ${CFDTOOLS_LIBRARY}")
        SET(CFDTOOLS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
MESSAGE("----fiun--> ${CFDTOOLS_LIBRARY}------<")
message(${CFDTOOLS_INCLUDE_DIR})
MARK_AS_ADVANCED(CFDTOOLS_LIBRARY)
MARK_AS_ADVANCED(CFDTOOLS_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${CFDTOOLS_INCLUDE_DIR})
