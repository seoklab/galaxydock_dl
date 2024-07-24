cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
project(GalaxyDistribution LANGUAGES C Fortran)

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(FATAL_ERROR "Only Intel Fortran is supported")
endif()

include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  "${PROJECT_SOURCE_DIR}/${CMAKE_INSTALL_BINDIR}")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS OFF)

string(APPEND CMAKE_C_FLAGS " -Wall -Wextra -pedantic")
string(APPEND CMAKE_Fortran_FLAGS " -fpp -heap-arrays -free \
-warn all -implicit-none -diag-enable=remark \
-diag-disable=7712,8889,10120,11000,11001,11006,11074,11076")
add_link_options(
  -static-intel
  -static-libgcc
  -static-libstdc++
)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O2 -no-prec-div -fp-model fast=2)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

if(NOT CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")
  add_link_options(-Xlinker --strip-all)
endif()

set(BLA_STATIC ON)
set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK REQUIRED)

file(GLOB_RECURSE GALAXY_SRCS "src/*.f90" "src/*.c")
list(FILTER GALAXY_SRCS EXCLUDE REGEX "mpi_[^/]*$")

add_executable(ligdock ${GALAXY_SRCS})
target_link_libraries(ligdock PRIVATE LAPACK::LAPACK)
target_include_directories(ligdock SYSTEM PRIVATE "$ENV{MKLROOT}/include/fftw")
