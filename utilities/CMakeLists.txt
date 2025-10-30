file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")


set(BOOST_INCLUDE_DIR ${DEFAULT_BOOST_ROOT}/include)
set(BOOST_LIBRARY_DIR ${DEFAULT_BOOST_ROOT}/stage/lib)


IF(THIRDPARTY_BUILD_VTK)

set(VTK_LIBRARIES
  ${VTK_LIBRARY}/libvtksys-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkCommonCore-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkCommonDataModel-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOLegacy-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOXML-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkParallelCore-9.4${LIB_EXT}
  ${VTK_LIBRARY}/libvtkIOParallelXML-9.4${LIB_EXT}
)

ELSE()
find_package(VTK)
target_link_libraries(${TEST_NAME} VTK::CommonCore VTK::CommonDataModel VTK::IOLegacy VTK::IOXML VTK::ParallelCore VTK::IOParallelXML)

ENDIF()


include_directories(${BOOST_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIR})
include_directories(${METIS_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${MMG_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${ESP_INCLUDE_DIR})

find_package(LAPACK)

add_subdirectory(useESP)
