set(TEST_NAME "tesselateRemainder")

add_executable(${TEST_NAME} main.cpp ${SRC_FILES})
target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

#INCLUDE (../../cmake/ThirdPartyOCCT.cmake)
#INCLUDE (../../cmake/ThirdPartyTetGen.cmake)


include_directories(${OCCT_INCLUDE_DIR})
include_directories(${TETGEN_INCLUDE_DIR})
#target_link_libraries(${TEST_NAME}
#    TKernel
#    TKMath
#    TKBRep
#    TKGeomBase
#    TKGeomAlgo
#    TKTopAlgo
#    TKPrim
#    TKBO
#    TKMesh
#    TKSTEP
#    TKSTEPBase
#    TKSTEPAttr
#    TKSTEP209
#    TKSTL
#    TKXSBase
#)

target_link_libraries(${TEST_NAME}
	${OCCT_LIBRARY}/libTKG2d${LIB_EXT}
	${OCCT_LIBRARY}/libTKG3d${LIB_EXT}
	${OCCT_LIBRARY}/libTKShHealing${LIB_EXT}
	${OCCT_LIBRARY}/libTKernel${LIB_EXT}
	${OCCT_LIBRARY}/libTKMath${LIB_EXT}
	${OCCT_LIBRARY}/libTKBRep${LIB_EXT}
	${OCCT_LIBRARY}/libTKGeomBase${LIB_EXT}
	${OCCT_LIBRARY}/libTKGeomAlgo${LIB_EXT}
	${OCCT_LIBRARY}/libTKTopAlgo${LIB_EXT}
	${OCCT_LIBRARY}/libTKBO${LIB_EXT}
	${OCCT_LIBRARY}/libTKMesh${LIB_EXT}
	${OCCT_LIBRARY}/libTKDESTEP${LIB_EXT}
	${OCCT_LIBRARY}/libTKDEIGES${LIB_EXT}
	${OCCT_LIBRARY}/libTKDESTL${LIB_EXT}
	${OCCT_LIBRARY}/libTKXSBase${LIB_EXT}
	${OCCT_LIBRARY}/libTKPrim${LIB_EXT}
)
target_link_libraries(${TEST_NAME} ${ESP_LIBRARY}/libegads${LIB_EXT})
target_link_libraries(${TEST_NAME} ${TETGEN_LIBRARY})
target_link_libraries(${TEST_NAME} ${LAPACK_LIBRARIES})
target_link_libraries(${TEST_NAME} ${METIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${PARMETIS_LIBRARY})
target_link_libraries(${TEST_NAME} ${HDF5_LIBRARY})
target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})

target_link_libraries(${TEST_NAME} ${MMG_LIBRARY}/libmmg3d${LIB_EXT})
target_link_libraries(${TEST_NAME} ${PARMMG_LIBRARY}/libparmmg${LIB_EXT})

#add_test(NAME tesselateRemainder COMMAND tesselateRemainder)


