# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build

# Utility rule file for stencil-classes.cc.

# Include any custom commands dependencies for this target.
include blitz/array/CMakeFiles/stencil-classes.cc.dir/compiler_depend.make

# Include the progress variables for this target.
include blitz/array/CMakeFiles/stencil-classes.cc.dir/progress.make

blitz/array/CMakeFiles/stencil-classes.cc:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array && /home/johannes/.pyenv/shims/python3.9 /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/../generate/genstencils.py stencil-classes.cc MAIN_DEPENDENCY /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array/../generate/genstencils.py

stencil-classes.cc: blitz/array/CMakeFiles/stencil-classes.cc
stencil-classes.cc: blitz/array/CMakeFiles/stencil-classes.cc.dir/build.make
.PHONY : stencil-classes.cc

# Rule to build all files generated by this target.
blitz/array/CMakeFiles/stencil-classes.cc.dir/build: stencil-classes.cc
.PHONY : blitz/array/CMakeFiles/stencil-classes.cc.dir/build

blitz/array/CMakeFiles/stencil-classes.cc.dir/clean:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array && $(CMAKE_COMMAND) -P CMakeFiles/stencil-classes.cc.dir/cmake_clean.cmake
.PHONY : blitz/array/CMakeFiles/stencil-classes.cc.dir/clean

blitz/array/CMakeFiles/stencil-classes.cc.dir/depend:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2 /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/array /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/array/CMakeFiles/stencil-classes.cc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : blitz/array/CMakeFiles/stencil-classes.cc.dir/depend

