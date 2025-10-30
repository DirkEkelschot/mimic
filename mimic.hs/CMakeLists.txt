
#file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

set(TEST_NAME "mimic.hs")

add_executable(${TEST_NAME} main.cpp ${SRC_FILES})
target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR})


target_link_libraries(${TEST_NAME} ${LAPACK_LIBRARIES})
target_link_libraries(${TEST_NAME} ${METIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${PARMETIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${HDF5_LIBRARY})
target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})
target_link_libraries(${TEST_NAME} ${MMG_LIBRARY}/libmmg3d${LIB_EXT})
target_link_libraries(${TEST_NAME} ${PARMMG_LIBRARY}/libparmmg${LIB_EXT})

add_test(NAME mimic.su2 COMMAND mimic.hs)


