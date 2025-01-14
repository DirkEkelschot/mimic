# Configuration file for MacOSX

cmake_minimum_required(VERSION 3.1...3.16)

#set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_STANDARD 17)
INCLUDE(CMakeDependentOption)
#set(CMAKE_C_COMPILER icpc)

# set the project name and version
project(MIMIC VERSION 1.0)
# specify the C++ standard
include(./NAS_Dirk_work.cmake)
set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")
file(GLOB SRC CONFIGURE_DEPENDS "src/*.cpp" "main.cpp")
find_package(MPI REQUIRED)
#find_package(Boost COMPONENTS system filesystem REQUIRED)
#set(CMAKE_CXX_FLAGS "-Wl,-no_compact_unwind")

add_executable(computeMetricDistance ${SRC})

find_package(LAPACK)
#find_package(HDF5)
INCLUDE (cmake/adaptCommon.cmake)
INCLUDE (cmake/ThirdPartyGSL.cmake)
#find_package(Boost COMPONENTS system filesystem REQUIRED)
set(Boost_INCLUDE_DIR /nobackupp19/dekelsch/boost_1_71_0/boost-install/include)
set(Boost_LIBRARY_DIR /nobackupp19/dekelsch/boost_1_71_0/boost-install/stage/lib)

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS})


link_directories(${Boost_LIBRARY_DIR})
target_link_libraries(computeMetricDistance ${GSL_LIBRARY})
target_link_libraries(computeMetricDistance ${MPI_CXX_LIBRARIES})
target_link_libraries(computeMetricDistance ${LAPACK_LIBRARIES})
#target_link_libraries(partTetra ${Boost_LIBRARIES})
#target_link_libraries(partTetra ${Boost_LIBRARY_DIR})
#====== Adding the external libraries here ======



