SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

#FIND_LIBRARY(EIGEN_LIBRARY NAMES eigen)
#FIND_PATH(EIGEN_INCLUDE_DIR EIGEN.h)
IF(EIGEN_LIBRARY AND EIGEN_INCLUDE_DIR)
        SET(BUILD_EIGEN OFF)
ELSE()
        SET(BUILD_EIGEN ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_EIGEN
        "Build EIGEN library from ThirdParty." ${BUILD_EIGEN})

IF (THIRDPARTY_BUILD_EIGEN)
        INCLUDE(ExternalProject)
        
        EXTERNALPROJECT_ADD(
                eigen-3.4.0
                URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/eigen-3.4.0
                BINARY_DIR "" # Not needed for header-only libraries
                TMP_DIR ${TPBUILD}/eigen-3.4.0-tmp
                INSTALL_DIR ${TPDIST}
                CONFIGURE_COMMAND "" # No configure step needed
                BUILD_COMMAND "" # No build step needed
                INSTALL_COMMAND "" # No install step needed
                )
        SET(EIGEN_INCLUDE_DIR ${TPSRC}/eigen-3.4.0 CACHE FILEPATH
            "EIGEN include" FORCE)
        ADD_DEFINITIONS(-DEIGEN_HAS_DEINITIALIZE)
        message(${EIGEN_INCLUDE_DIR})
        MARK_AS_ADVANCED(EIGEN_INCLUDE_DIR)
        INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})
ENDIF()
