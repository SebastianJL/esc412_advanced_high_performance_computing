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

# Include any dependencies generated for this target.
include blitz/generate/CMakeFiles/genvecbfn.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include blitz/generate/CMakeFiles/genvecbfn.dir/compiler_depend.make

# Include the progress variables for this target.
include blitz/generate/CMakeFiles/genvecbfn.dir/progress.make

# Include the compile flags for this target's objects.
include blitz/generate/CMakeFiles/genvecbfn.dir/flags.make

blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o: blitz/generate/CMakeFiles/genvecbfn.dir/flags.make
blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o: /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genvecbfn.cpp
blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o: blitz/generate/CMakeFiles/genvecbfn.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o -MF CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o.d -o CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o -c /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genvecbfn.cpp

blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/genvecbfn.dir/genvecbfn.cpp.i"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genvecbfn.cpp > CMakeFiles/genvecbfn.dir/genvecbfn.cpp.i

blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/genvecbfn.dir/genvecbfn.cpp.s"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genvecbfn.cpp -o CMakeFiles/genvecbfn.dir/genvecbfn.cpp.s

# Object files for target genvecbfn
genvecbfn_OBJECTS = \
"CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o"

# External object files for target genvecbfn
genvecbfn_EXTERNAL_OBJECTS =

blitz/generate/genvecbfn: blitz/generate/CMakeFiles/genvecbfn.dir/genvecbfn.cpp.o
blitz/generate/genvecbfn: blitz/generate/CMakeFiles/genvecbfn.dir/build.make
blitz/generate/genvecbfn: blitz/generate/CMakeFiles/genvecbfn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable genvecbfn"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/genvecbfn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
blitz/generate/CMakeFiles/genvecbfn.dir/build: blitz/generate/genvecbfn
.PHONY : blitz/generate/CMakeFiles/genvecbfn.dir/build

blitz/generate/CMakeFiles/genvecbfn.dir/clean:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && $(CMAKE_COMMAND) -P CMakeFiles/genvecbfn.dir/cmake_clean.cmake
.PHONY : blitz/generate/CMakeFiles/genvecbfn.dir/clean

blitz/generate/CMakeFiles/genvecbfn.dir/depend:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2 /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate/CMakeFiles/genvecbfn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : blitz/generate/CMakeFiles/genvecbfn.dir/depend

