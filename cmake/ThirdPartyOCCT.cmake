#FIND_LIBRARY(VTK_LIBRARY NAMES VTK)
#FIND_PATH(TETGEN_INCLUDE_DIR VTK.h)

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)



IF(OCCT_LIBRARY AND OCCT_INCLUDE_DIR)
        SET(BUILD_OCCT OFF)
ELSE()
        SET(BUILD_OCCT ON)
ENDIF()



OPTION(THIRDPARTY_BUILD_OpenCascade
        "Build OpenCascade library from ThirdParty." ${BUILD_OCCT})


MESSAGE("BUILDING OCCT THIRDPARTY_BUILD_OpenCascade " ${BUILD_OCCT} ${THIRDPARTY_BUILD_OpenCascade})


IF (THIRDPARTY_BUILD_OpenCascade)
        INCLUDE(ExternalProject)
        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
                MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)
        
        EXTERNALPROJECT_ADD(
                OCCT-7_9_1
                PREFIX ${CMAKE_SOURCE_DIR}/Thirdparty
                URL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_9_1.tar.gz
		#URL_MD5 90d070d56d2607e200df0b56fde92cd4
		STAMP_DIR ${TPBUILD}/stamp
            	DOWNLOAD_DIR ${TPSRC}
            	SOURCE_DIR ${TPSRC}/OCCT-7_9_1
            	BINARY_DIR ${TPBUILD}/OCCT-7_9_1
            	TMP_DIR ${TPBUILD}/OCCT-7_9_1-tmp
            	INSTALL_DIR ${TPDIST}
		CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER=${MPI_C} -DCMAKE_CXX_COMPILER=${MPI_CXX} 
                -DCMAKE_CXX_COMPILER:FILEPATH=${MPI_CXX}
		-DCMAKE_C_COMPILER:FILEPATH=${MPI_C}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DVTK_ENABLE_PARALLEL=ON
		-DZLIB_DIR=${TPDIST}
                -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
                -DCMAKE_INSTALL_LIBDIR=lib
                # CRITICAL: Enable all required modules
                -DBUILD_MODULE_FoundationClasses=ON
                -DBUILD_MODULE_ModelingData=ON
                -DBUILD_MODULE_ModelingAlgorithms=ON
                -DBUILD_MODULE_Visualization=ON
                -DBUILD_MODULE_ApplicationFramework=ON
                -DBUILD_MODULE_DataExchange=ON          # THIS IS THE KEY ONE!
                # Optional but recommended
                -DBUILD_MODULE_Draw=OFF                  # Drawing/testing (not needed)
                -DBUILD_DOC_Overview=OFF                 # Documentation (not needed)
		${TPSRC}/OCCT-7_9_1
                BUILD_COMMAND make -j14
                DOWNLOAD_EXTRACT_TIMESTAMP TRUE
                )
        SET(OCCT_INCLUDE_DIR ${TPDIST}/include/opencascade CACHE FILEPATH
            "OCCT include directory" FORCE)
        SET(OCCT_LIBRARY ${TPDIST}/lib CACHE FILEPATH
            "OCCT lib directory" FORCE)
	SET(OCCT_CONFIG_INCLUDE_DIR ${TPINC})

ENDIF()
#SET(VTK_LIB ${TPDIST}/lib)

MESSAGE("OCCT_INCLUDE_DIR OCCT_INCLUDE_DIROCCT_INCLUDE_DIROCCT_INCLUDE_DIR.   " ${OCCT_INCLUDE_DIR})

#ADD_DEPENDENCIES(mimic OCCT-7_9_1)
#SET(VTK_DIR   ${TPDIST}lib/cmake/vtk-9.4)
#SET(VTK_LIBRARIES ${TPDIST}/lib CACHE FILEPATH
#            "VTK lib directory" FORCE)
MARK_AS_ADVANCED(OCCT_LIBRARY)
INCLUDE_DIRECTORIES(${OCCT_INCLUDE_DIR})
MARK_AS_ADVANCED(OCCT_INCLUDE_DIR)



