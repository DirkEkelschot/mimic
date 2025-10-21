#FIND_LIBRARY(PARMETIS_LIBRARY NAMES parmetis)
#FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(PARMETIS_LIBRARY AND PARMETIS_INCLUDE_DIR)
        SET(BUILD_PARMETIS OFF)
ELSE()
        SET(BUILD_PARMETIS ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_PARMETIS
        "Build TetGen library from ThirdParty." ${BUILD_PARMETIS})
message("CMAKE_BINARY_DIR=",${CMAKE_BINARY_DIR})
IF (THIRDPARTY_BUILD_PARMETIS)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build parmetis.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                parmetis
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty                
		GIT_REPOSITORY https://github.com/KarypisLab/ParMETIS.git
                GIT_TAG main
                STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/parmetis
            	BINARY_DIR ${TPBUILD}/parmetis 
            	TMP_DIR ${TPBUILD}/parmetis-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DGKLIB_PATH=${TPDIST}/gklib
                -DMETIS_PATH=${TPDIST}/metis
		${TPSRC}/parmetis
		)
        THIRDPARTY_LIBRARY(PARMETIS_LIBRARY STATIC parmetis
            DESCRIPTION "PARMETIS library")

        MESSAGE(STATUS "Build PARMETIS: ${PARMETIS_LIBRARY}")
        SET(PARMETIS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()
#ADD_DEPENDENCIES(thirdparty PARMETIS)
MARK_AS_ADVANCED(PARMETIS_LIBRARY)
MARK_AS_ADVANCED(PARMETIS_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})


