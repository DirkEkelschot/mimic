SET(TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
SET(TPBUILD ${CMAKE_BINARY_DIR}/ThirdParty)
SET(TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
SET(GKLIB_VERSION master) # Set as needed



find_library(GKLIB_LIBRARY libGKlib.a ${TPDIST}/gklib/lib)

#if(GKLIB_LIBRARY)
#    set(BUILD_GKLIB OFF)
#else()
#    set(BUILD_GKLIB ON)
#endif()

set(BUILD_GKLIB ON)
IF (BUILD_GKLIB)
    INCLUDE(ExternalProject)
    UNSET(MAKE_EXE CACHE)
    FIND_PROGRAM(MAKE_EXE NAMES gmake make)
    IF(NOT MAKE_EXE)
        MESSAGE(FATAL_ERROR
        "'make' tool for building GKlib not found. Cannot build GKlib.")
    ENDIF()
    MARK_AS_ADVANCED(MAKE_EXE)

    SET(GKLIB_BUILD_DIR ${TPBUILD}/GKlib)
    SET(GKLIB_INSTALL_DIR ${TPDIST}/gklib)

    ExternalProject_Add(
        gklib
        GIT_REPOSITORY https://github.com/KarypisLab/GKlib.git
        GIT_TAG ${GKLIB_VERSION}
        SOURCE_DIR ${TPSRC}/GKlib
        BINARY_DIR ${TPSRC}/GKlib
        INSTALL_DIR ${TPDIST}/gklib
        CONFIGURE_COMMAND make config prefix=${TPDIST}/gklib
        BUILD_COMMAND make
        INSTALL_COMMAND make install prefix=${TPDIST}/gklib
        && ${CMAKE_COMMAND} -E copy ${TPDIST}/gklib/lib64/libGKlib.a ${TPDIST}/gklib/lib/libGKlib.a
    )

    # Mark installed library and include as external variables
    SET(GKLIB_LIBRARY ${GKLIB_INSTALL_DIR}/lib/libGKlib.a)
    SET(GKLIB_INCLUDE_DIR ${GKLIB_INSTALL_DIR}/include)
    INCLUDE_DIRECTORIES(${GKLIB_INCLUDE_DIR})
    MARK_AS_ADVANCED(GKLIB_LIBRARY GKLIB_INCLUDE_DIR)
    MESSAGE(STATUS "Build GKlib: ${GKLIB_LIBRARY}")
    set(GKLIB_DEP gklib)
ELSE()
    set(GKLIB_DEP)
ENDIF()

