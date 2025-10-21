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
MESSAGE("TPSRC" ${TPSRC}) 


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

        EXTERNALPROJECT_ADD(
                metis
                GIT_REPOSITORY https://github.com/KarypisLab/METIS.git
                GIT_TAG master
                SOURCE_DIR ${TPSRC}/metis
                BINARY_DIR ${TPSRC}/metis
            	INSTALL_DIR ${TPDIST}
                CONFIGURE_COMMAND make config shared=1 prefix=${TPDIST}/metis gklib_path=${TPDIST}/gklib
                BUILD_COMMAND make CFLAGS="-I${TPDIST}/gklib/include" -C ${TPSRC}/metis gklib_path=${TPDIST}/gklib
                INSTALL_COMMAND make CFLAGS="-I${TPDIST}/gklib/include" -C ${TPSRC}/metis install gklib_path=${TPDIST}/gklib
                DEPENDS gklib
		)
        #THIRDPARTY_LIBRARY(METIS_LIBRARY STATIC metis
        #    DESCRIPTION " library")

        #MESSAGE(STATUS "Build METIS: ${METIS_LIBRARY}")
        #SET(METIS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
#ADD_DEPENDENCIES(thirdparty METIS)
#MARK_AS_ADVANCED(METIS_LIBRARY)
#MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
#INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})





