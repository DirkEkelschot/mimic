#FIND_LIBRARY(MMG_LIBRARY NAMES mmg)
#FIND_PATH(MMG_INCLUDE_DIR mmg.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(MMG_LIBRARY AND MMG_INCLUDE_DIR)
        SET(BUILD_MMG OFF)
ELSE()
        SET(BUILD_MMG ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_MMG
        "Build TetGen library from ThirdParty." ${BUILD_MMG})
message("CMAKE_BINARY_DIR=",${CMAKE_BINARY_DIR})
IF (THIRDPARTY_BUILD_MMG)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build mmg.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                mmg
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                #GIT_REPOSITORY https://github.com/MmgTools/mmg.git
                #GIT_TAG master
		GIT_REPOSITORY https://github.com/MmgTools/mmg.git
    		GIT_TAG v5.8.0
                STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/mmg
            	BINARY_DIR ${TPBUILD}/mmg 
            	TMP_DIR ${TPBUILD}/mmg-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
		-DMMG_INSTALL_PRIVATE_HEADERS=ON
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
		${TPSRC}/mmg
		)
        THIRDPARTY_LIBRARY(MMG_LIBRARY STATIC mmg
            DESCRIPTION "MMGS library")
	
	SET(MMG_DIR ${TPDIST} CACHE FILEPATH
            "MMG DIR directory" FORCE)
        SET(MMG_INCLUDE_DIRS ${TPDIST}/include CACHE FILEPATH
            "MMG include directory" FORCE)
	SET(MMG_INCLUDE_LIB ${TPDIST}/lib CACHE FILEPATH
            "MMG lib directory" FORCE)
        MESSAGE(STATUS "Build MMG: ${MMG_LIBRARY}")
        SET(MMG_CONFIG_INCLUDE_DIR ${TPINC})
	SET(MMG_INCLUDE_DIRS ${TPDIST}/include)
	SET(MMGS_INCLUDE_DIR ${TPDIST}/include/mmg/mmgs)
	SET(MMG3D_INCLUDE_DIR ${TPDIST}/include/mmg/mmg3d)
ENDIF()
message("MMG_DIR======== ",${MMG_DIR})
#ADD_DEPENDENCIES(thirdparty MMG)
MARK_AS_ADVANCED(MMG_LIBRARY)
MARK_AS_ADVANCED(MMG_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIRS})


