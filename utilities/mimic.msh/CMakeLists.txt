
file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

set(TEST_NAME "mimic.msh")

add_executable(${TEST_NAME} main.cpp ${SRC_FILES})
target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

set(BOOST_INCLUDE_DIR ${DEFAULT_BOOST_ROOT}/include)
set(BOOST_LIBRARY_DIR ${DEFAULT_BOOST_ROOT}/stage/lib)


IF(NOT VTK_DIR)

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

MESSAGE("VTK_LIBRARIES " ${PARMETIS_INCLUDE_DIR})

include_directories(${BOOST_INCLUDE_DIR})
include_directories(${PARMETIS_INCLUDE_DIR})
include_directories(${METIS_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${MMG_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})

find_package(LAPACK)
target_link_libraries(${TEST_NAME} ${LAPACK_LIBRARIES})
target_link_libraries(${TEST_NAME} ${METIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${PARMETIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${HDF5_LIBRARY})
target_link_libraries(${TEST_NAME} ${MMG_LIBRARY}/libmmg3d${LIB_EXT})
target_link_libraries(${TEST_NAME} ${PARMMG_LIBRARY}/libparmmg${LIB_EXT})

add_test(NAME mimic.msh COMMAND mimic.msh)


