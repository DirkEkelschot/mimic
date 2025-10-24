

SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)

IF(ESP_LIBRARY AND ESP_INCLUDE_DIR)
        SET(BUILD_ESP OFF)
ELSE()
        SET(BUILD_ESP ON)
ENDIF()

OPTION(THIRDPARTY_BUILD_ESP
        "Build TetGen library from ThirdParty." ${BUILD_ESP})

IF (THIRDPARTY_BUILD_ESP)
        INCLUDE(ExternalProject)
#        UNSET(PATCH CACHE)
#        FIND_PROGRAM(PATCH patch)
#        IF(NOT PATCH)
#                MESSAGE(FATAL_ERROR
#                "'patch' tool for modifying files not found. Cannot build ESP.")
#        ENDIF()
#        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
                ESP
		URL https://acdl.mit.edu/ESP/ESP.tgz
                STAMP_DIR ${TPBUILD}/stamp
                DOWNLOAD_DIR ${TPSRC}
                SOURCE_DIR ${TPSRC}/ESP
                #BINARY_DIR ${TPBUILD}/ESP
                TMP_DIR ${TPBUILD}/ESP-tmp
                BUILD_IN_SOURCE 1
		CONFIGURE_COMMAND cd ${TPSRC}/ESP/config && ./makeEnv ${TPDIST}
		BUILD_COMMAND bash -c "cd ${TPSRC}/ESP && source ESPenv.sh && cd src && make && cp -r ${TPSRC}/ESP/include/*.h ${TPDIST}/include && cp -r ${TPSRC}/ESP/lib/* ${TPDIST}/lib"
                INSTALL_COMMAND cmake -E echo "Skipping install step for ESP"
#		INSTALL_COMMAND make install
                ${TPSRC}/ESP
                )
        THIRDPARTY_LIBRARY(ESP_LIBRARY STATIC ESP
            DESCRIPTION "ESP library")
	
	SET(ESP_DIR ${TPDIST} CACHE FILEPATH
            "ESP DIR directory" FORCE)
        SET(ESP_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "ESP include directory" FORCE)
	SET(ESP_LIBRARY ${TPDIST}/lib CACHE FILEPATH
            "ESP lib directory" FORCE)
        MESSAGE(STATUS "Build ESP: ${ESP_LIBRARY}")
        SET(ESP_CONFIG_INCLUDE_DIR ${TPINC})
ENDIF()

MESSAGE("ESP_LIBRARY   " ${ESP_LIBRARY})
ADD_DEPENDENCIES(mimic ESP)
MARK_AS_ADVANCED(ESP_LIBRARY)
MARK_AS_ADVANCED(ESP_INCLUDE_DIR)
INCLUDE_DIRECTORIES(${ESP_INCLUDE_DIR})


