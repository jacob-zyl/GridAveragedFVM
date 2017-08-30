cmake_minimum_required(VERSION 3.5)

project (HighOrderScheme)
enable_language (Fortran)

set (CMAKE_BUILD_TYPE DEBUG)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack")
set (CMAKE_Fortran_COMPILER_ID "GNU")

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set (CMAKE_Fortran_COMPILER "gfortran")
  set (CMAKE_C_COMPILER "gcc")
  set (CMAKE_CXX_COMPILER "g++")
  set (dialect "-ffree-form -std=f2008 -fimplicit-none")
  set (bounds "-fbounds-check")
  set (warning "-Wall -Wextra")
  set (precision "-fdefault-real-8")
  set (performance "-Ofast")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set (CMAKE_Fortran_COMPILER "ifort")
  set (CMAKE_C_COMPILER "icc")
  set (CMAKE_CXX_COMPILER "icpc")
  set (dialect "-stand f08 -free -implicitnone")
  set (bounds "-check bounds")
  set (warning "-warn")
  set (precision "-r8")
  set (performance "-Ofast")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set (CMAKE_Fortran_COMPILER "pgfortran")
  set (CMAKE_C_COMPILER "pgcc")
  set (CMAKE_CXX_COMPILER "pgc++")
  set (dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
  set (bounds "-C")
  set (warning "")
  set (precision "-r8")
  set (performance "-fast -Mipa=fast,inline")
endif ()

set (CMAKE_Fortran_FLAGS_DEBUG
  "${CMAKE_Fortran_FLAGS} ${dialect} ${precision} ${bounds} ${warning}")
set (CMAKE_Fortran_FLAGS_RELAESE
  "${CMAKE_Fortran_FLAGS} ${dialect} ${precision} ${performance}")

#
# Compile.
#
file (GLOB_RECURSE sources  src/main.f90 src/scheme.f90 src/constant.f90)
add_executable (prog ${sources})
