cmake_minimum_required(VERSION 2.8)

project(FLUSH)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -Wall -fno-strict-aliasing -mssse3" )
#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall  -Wshadow -Wno-error=shadow -Wno-error=unused-function -fomit-frame-pointer -ffast-math -fvisibility=hidden -fvisibility-inlines-hidden -fPIC  -Werror " )

message(${CMAKE_HOST_SYSTEM_VERSION})

if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^4.4" )
  message( "Hopefully, a Ubuntu 16.04 machine or similar, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine with an ancient compiler." )
endif()

#############################################################################
## required packages
#############################################################################

find_package(Boost COMPONENTS program_options regex date_time system filesystem REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

#############################################################################

set(EXECUTABLE_OUTPUT_PATH ${FLUSH_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

set(DACLIB_SRCS2 build_time.cc get_cwd.cc)

set(DACLIB_SRCS3 build_time.cc
get_cwd.cc
mpi_string_subs.cc)

set(FP_SRCS FingerprintBase.cc
HashedFingerprint.cc
NotHashedFingerprint.cc)

set(DACLIB_INCS3
ByteSwapper.H
FileExceptions.H
FingerprintBase.H
HashedFingerprint.H
MagicInts.H
NotHashedFingerprint.H)

set(FP_INCS FingerprintBase.H
HashedFingerprint.H
NotHashedFingerprint.H)

#############################################################################
## satan, cluster, amtec, subset_fp_file, merge_fp_files, cad, histogram
#############################################################################

add_executable(satan satan.cc
SatanSettings.cc
${FP_SRCS} ${DACLIB_SRCS3} ${DACLIB_INCS3} ${FP_INCS})

target_link_libraries(satan ${LIBS} ${Boost_LIBRARIES}
${MPI_LIBRARIES} z)

add_executable(cluster cluster.cc
ClusterSettings.cc
${FP_SRCS} ${DACLIB_SRCS3})

target_link_libraries(cluster ${LIBS} ${Boost_LIBRARIES}
  ${MPI_LIBRARIES} z)

add_executable(amtec amtec.cc
AmtecSettings.cc
${FP_SRCS} build_time.cc)

target_link_libraries(amtec ${LIBS} ${Boost_LIBRARIES} z)

add_executable(subset_fp_file subset_fp_file.cc
${FP_SRCS} build_time.cc)

target_link_libraries(subset_fp_file ${LIBS} ${Boost_LIBRARIES} z)

add_executable(reverse_fp_file reverse_fp_file.cc
${FP_SRCS} build_time.cc)

target_link_libraries(reverse_fp_file ${LIBS} ${Boost_LIBRARIES} z)

add_executable(merge_fp_files merge_fp_files.cc
${FP_SRCS} build_time.cc)

target_link_libraries(merge_fp_files ${LIBS} ${Boost_LIBRARIES} z)

add_executable(cad cad.cc
CadSettings.cc ${DACLIB_SRCS2}
${FP_SRCS} build_time.cc)

target_link_libraries(cad ${LIBS} ${Boost_LIBRARIES} z)

add_executable(histogram histogram.cc
${FP_SRCS} ${DACLIB_SRCS2})
target_link_libraries(histogram ${LIBS} ${Boost_LIBRARIES} z)
