#FIND_LIBRARY(METIS_LIBRARY NAMES metis)
#FIND_PATH(METIS_INCLUDE_DIR metis.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(METIS_LIBRARY AND METIS_INCLUDE_DIR)
        SET(BUILD_METIS OFF)
ELSE()
        SET(BUILD_METIS ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_METIS
        "Build TetGen library from ThirdParty." ${BUILD_METIS})


find_library(GKLIB_LIBRARY libGKlib.a ${TPDIST}/gklib/lib)

MESSAGE("TPBUILD" ${TPBUILD}) 
#MESSAGE("TPSRC --- CONFIGURE_COMMAND " ${TPSRC} " CONFIGURE_COMMAND " ${CONFIGURE_COMMAND}) 
if (CMAKE_VERSION VERSION_GREATER "3.23.4")
    cmake_policy(SET CMP0135 NEW)
endif()

IF (THIRDPARTY_BUILD_METIS)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build METIS.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)
        
        SET(METIS_BUILD_DIR ${TPBUILD}/metis)
        SET(METIS_INSTALL_DIR ${TPDIST}/metis)

        ExternalProject_Add(
                        metis
                        URL https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz
                        URL_MD5 5465e67079419a69e0116de24fce58fe
                        UPDATE_COMMAND ""
                        CONFIGURE_COMMAND ${CMAKE_MAKE_PROGRAM} config prefix=../../../Metis-install
                        BUILD_IN_SOURCE 1
                        BUILD_COMMAND unset MFLAGS && unset MAKEFLAGS && unset MAKELEVEL && ${CMAKE_MAKE_PROGRAM}
                        INSTALL_COMMAND unset MFLAGS && unset MAKEFLAGS && unset MAKELEVEL && ${CMAKE_MAKE_PROGRAM} install
                        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
                        )
        THIRDPARTY_LIBRARY(METIS_LIBRARY STATIC metis
            DESCRIPTION " library")

        MESSAGE(STATUS "Build METIS: ${METIS_LIBRARY}")
        SET(METIS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
#ADD_DEPENDENCIES(thirdparty METIS)
MARK_AS_ADVANCED(METIS_LIBRARY)
MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})





