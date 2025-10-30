#FIND_LIBRARY(VTK_LIBRARY NAMES VTK)
#FIND_PATH(TETGEN_INCLUDE_DIR VTK.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)



IF(VTK_LIBRARY AND VTK_INCLUDE_DIR)
        SET(BUILD_VTK OFF)
ELSE()
        SET(BUILD_VTK ON)
ENDIF()



OPTION(THIRDPARTY_BUILD_VTK
        "Build TetGen library from ThirdParty." ${BUILD_VTK})

IF (THIRDPARTY_BUILD_VTK)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                VTK-9.4.2
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL https://vtk.org/files/release/9.4/VTK-9.4.2.tar.gz
		URL_MD5 90d070d56d2607e200df0b56fde92cd4
		STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/VTK-9.4.2
            	BINARY_DIR ${TPBUILD}/VTK-9.4.2
            	TMP_DIR ${TPBUILD}/VTK-9.4.2-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DVTK_ENABLE_PARALLEL=ON
		-DZLIB_DIR=${TPDIST}
                -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
                -DCMAKE_INSTALL_LIBDIR=lib
                ${TPSRC}/VTK-9.4.2
                BUILD_COMMAND make -j14
                DOWNLOAD_EXTRACT_TIMESTAMP TRUE
                )
        SET(VTK_INCLUDE_DIRS ${TPDIST}/include/vtk-9.4 CACHE FILEPATH
            "VTK include directory" FORCE)
        SET(VTK_LIBRARY ${TPDIST}/lib CACHE FILEPATH
            "VTK lib directory" FORCE)
	SET(VTK_CONFIG_INCLUDE_DIR ${TPINC})

ENDIF()
#SET(VTK_LIB ${TPDIST}/lib)


#ADD_DEPENDENCIES(mimic VTK-9.4.2)
SET(VTK_DIR   ${TPDIST}lib/cmake/vtk-9.4)
SET(VTK_LIBRARIES ${TPDIST}/lib CACHE FILEPATH
            "VTK lib directory" FORCE)
MARK_AS_ADVANCED(VTK_LIBRARIES)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
MARK_AS_ADVANCED(VTK_INCLUDE_DIRS)



