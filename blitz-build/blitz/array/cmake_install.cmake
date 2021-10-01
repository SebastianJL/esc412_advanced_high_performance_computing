# Install script for directory: /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/johannes/uni/msc/uzh/HS21/ahpc/project")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/blitz/array" TYPE FILE FILES
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/asexpr.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/asexpr.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/cartesian.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/cgsolve.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/complex.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/convolve.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/convolve.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/cycle.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/domain.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/et.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/expr.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/expr.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/fastiter.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/funcs.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/functorExpr.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/geometry.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/indirect.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/interlace.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/io.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/iter.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/map.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/methods.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/misc.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/multi.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/newet-macros.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/newet.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/ops.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/ops.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/reduce.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/reduce.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/resize.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/shape.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/slice.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/slicing.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/stencil-et.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/stencil-et-macros.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/stencilops.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/stencils.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/stencils.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/storage.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/where.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/zip.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array/bops.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array/uops.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array/stencil-classes.cc"
    )
endif()

