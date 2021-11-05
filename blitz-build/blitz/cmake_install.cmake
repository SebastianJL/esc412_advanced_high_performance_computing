# Install script for directory: /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/blitz" TYPE FILE FILES
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array-impl.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/bench.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/bench.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/benchext.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/benchext.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/blitz.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/bounds.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/bzdebug.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/bzconfig.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/compiler.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/constpointerstack.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/etbase.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/et-forward.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/funcs.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/globeval.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/indexexpr.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/indexmap-forward.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/levicivita.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/limits-hack.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/listinit.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/memblock.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/memblock.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/minmax.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/numinquire.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/numtrait.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/ops.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/prettyprint.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/promote.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/range.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/range.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/ranks.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/reduce.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/shapecheck.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/simdtypes.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tau.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/timer.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinymat2.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinymat2.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinymat2io.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinyvec2.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinyvec2.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tinyvec2io.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tm2fastiter.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tmevaluate.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tv2fastiter.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tvevaluate.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/traversal.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/traversal.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tuning.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tvcross.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/tvecglobs.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/update.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/wrap-climits.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/vecbops.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/vecuops.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/vecwhere.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/vecbfn.cc"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/matbops.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/matuops.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/mathfunc.h"
    "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/promote-old.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/blitz" TYPE FILE FILES "/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz//bzconfig.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate/cmake_install.cmake")
  include("/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/meta/cmake_install.cmake")
  include("/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array/cmake_install.cmake")

endif()

