#FIND_LIBRARY(PARMMG_LIBRARY NAMES parmmg)
#FIND_PATH(PARMMG_INCLUDE_DIR parmmg.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
IF(PARMMG_LIBRARY AND PARMMG_INCLUDE_DIR)
        SET(BUILD_PARMMG OFF)
ELSE()
        SET(BUILD_PARMMG ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_PARMMG
        "Build TetGen library from ThirdParty." ${BUILD_PARMMG})
message("MMG_DIR======== INHERERER ",${MMG_DIR})
IF (THIRDPARTY_BUILD_PARMMG)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                ParMmg
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                GIT_REPOSITORY https://github.com/MmgTools/ParMmg.git
                GIT_TAG master
                STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/ParMmg
            	BINARY_DIR ${TPBUILD}/ParMmg 
            	TMP_DIR ${TPBUILD}/ParMmg-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  
		-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
#		-DDOWNLOAD_MMG=ON
		-DDOWNLOAD_MMG=OFF -DMMG_DIR=${TPDIST}
                -DDOWNLOAD_METIS=OFF -DMETIS_DIR=${DEFAULT_METIS_ROOT}
                ${TPSRC}/ParMmg
                )
        THIRDPARTY_LIBRARY(PARMMG_LIBRARY STATIC ParMmg
            DESCRIPTION "ParMmg library")
        SET(PARMMG_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "ParMmg include directory" FORCE)
#        SET(MMG_INCLUDE_DIRS ${TPDIST}/ParMmg/Mmg-prefix/src/Mmg-build/include CACHE FILEPATH
#            "MMG include directory" FORCE)
#        SET(MMG3D_INCLUDE_DIR ${TPDIST}/ParMmg/Mmg-prefix/src/Mmg-build/include/mmg/mmg3d)
#	SET(MMG_INCLUDE_LIB ${TPDIST}/ParMmg/Mmg-prefix/src/Mmg-build/lib CACHE FILEPATH
#            "MMG lib directory" FORCE)
	MESSAGE(STATUS "Build ParMmg: ${PARMMG_LIBRARY}")
        SET(ParMMMG_CONFIG_INCLUDE_DIR ${TPINC})
	
ENDIF()
#ADD_DEPENDENCIES(thirdparty parmmg)
MARK_AS_ADVANCED(PARMMG_LIBRARY)
MARK_AS_ADVANCED(PARMMG_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${PARMMG_INCLUDE_DIR})
#MARK_AS_ADVANCED(MMG_LIBRARY)
#MARK_AS_ADVANCED(MMG_INCLUDE_DIR)
#INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIRS})

