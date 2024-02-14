MACRO(THIRDPARTY_LIBRARY varname)
    CMAKE_PARSE_ARGUMENTS(TPLIB "" "DESCRIPTION" "STATIC;SHARED" ${ARGN})

    IF(TPLIB_SHARED)
        IF(WIN32)
            # Ensure linking against .lib files on Windows
            SET(LIBTYPE "STATIC")
        ELSE()
            SET(LIBTYPE "SHARED")
        ENDIF()
        SET(TPLIBS ${TPLIB_SHARED})
    ELSEIF(TPLIB_STATIC)
        SET(LIBTYPE "STATIC")
        SET(TPLIBS ${TPLIB_STATIC})
    ENDIF()

    FOREACH (lib ${TPLIBS})
        LIST(APPEND tmplist "${TPDIST}/lib/${CMAKE_${LIBTYPE}_LIBRARY_PREFIX}${lib}${CMAKE_${LIBTYPE}_LIBRARY_SUFFIX}")
    ENDFOREACH()

    SET(${varname} ${tmplist} CACHE FILEPATH ${TPLIB_DESCRIPTION} FORCE)
    UNSET(tmplist)
    UNSET(LIBTYPE)
    UNSET(TPLIBS)
    UNSET(TPLIB_SHARED)
    UNSET(TPLIB_STATIC)
    UNSET(lib)
ENDMACRO()
