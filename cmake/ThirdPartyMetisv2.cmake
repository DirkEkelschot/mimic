SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
SET(GKLIB_VERSION master) # Set as needed


find_library(GKLIB_LIBRARY libGKlib.a ${TPDIST}/gklib/lib)

#if(METIS_LIBRARY)
#    set(BUILD_METIS OFF)
#else()
#    set(BUILD_METIS ON)
#endif()

set(BUILD_METIS ON)
IF (BUILD_METIS)
    INCLUDE(ExternalProject)
    UNSET(MAKE_EXE CACHE)
    FIND_PROGRAM(MAKE_EXE NAMES gmake make)
    IF(NOT MAKE_EXE)
        MESSAGE(FATAL_ERROR
        "'make' tool for building metis not found. Cannot build metis.")
    ENDIF()
    MARK_AS_ADVANCED(MAKE_EXE)

    SET(METIS_BUILD_DIR ${TPBUILD}/metis)
    SET(METIS_INSTALL_DIR ${TPDIST}/metis)

    ExternalProject_Add(
        metis
        URL https://karypis.github.io/glaros/files/sw/parmetis/parmetis-4.0.3.tar.gz
        URL_MD5 5465e67079419a69e0116de24fce58fe
        SOURCE_DIR ${TPSRC}/metis
        BINARY_DIR ${TPSRC}/metis
        INSTALL_DIR ${TPDIST}/metis
        CONFIGURE_COMMAND make config shared=1 prefix=${TPDIST}/metis
    )

    # Mark installed library and include as external variables
    SET(METIS_LIBRARY ${METIS_INSTALL_DIR}/lib/metis.a)
    SET(METIS_INCLUDE_DIR ${METIS_INSTALL_DIR}/include)
    INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
    MARK_AS_ADVANCED(METIS_LIBRARY METIS_INCLUDE_DIR)
    MESSAGE(STATUS "Build metis: ${METIS_LIBRARY}")
    set(METIS_DEP metis)
ELSE()
    set(METIS_DEP)
ENDIF()

