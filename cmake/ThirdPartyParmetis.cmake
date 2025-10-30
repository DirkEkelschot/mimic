

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(METIS_LIBRARY AND metis_INCLUDE_DIR)
        SET(BUILD_METIS OFF)
ELSE()
        SET(BUILD_METIS ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_metis
        "Build TetGen library from ThirdParty." ${BUILD_metis})

IF (THIRDPARTY_BUILD_METIS)
        INCLUDE(ExternalProject)
#        UNSET(PATCH CACHE)
#        FIND_PROGRAM(PATCH patch)
#        IF(NOT PATCH)
#                MESSAGE(FATAL_ERROR
#                "'patch' tool for modifying files not found. Cannot build metis.")
#        ENDIF()
#        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                metis
		GIT_REPOSITORY https://github.com/KarypisLab/METIS.git
                GIT_TAG master
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/metis
                #BINARY_DIR ${TPBUILD}/metis
                TMP_DIR ${TPBUILD}/metis-tmp
                BUILD_IN_SOURCE 1
		CONFIGURE_COMMAND make config cc=${MPICC} gklib_path=${TPDIST} prefix=${TPDIST}
		BUILD_COMMAND bash -c "make && make install"
                INSTALL_COMMAND cmake -E echo "Skipping install step for metis"
#		INSTALL_COMMAND make install
                ${TPSRC}/metis
                )
        THIRDPARTY_LIBRARY(metis_LIBRARY STATIC metis
            DESCRIPTION "metis library")
	
	SET(METIS_DIR ${TPDIST} CACHE FILEPATH
            "metis DIR directory" FORCE)
        SET(METIS_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "metis include directory" FORCE)
	SET(METIS_LIBRARY ${TPDIST}/lib CACHE FILEPATH
            "metis lib directory" FORCE)
        MESSAGE(STATUS "Build metis: ${metis_LIBRARY}")
        SET(METIS_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()

MESSAGE("metis_LIBRARY   " ${METIS_LIBRARY})
ADD_DEPENDENCIES(mimic METIS)
MARK_AS_ADVANCED(METIS_LIBRARY)
MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${metis_INCLUDE_DIR})


