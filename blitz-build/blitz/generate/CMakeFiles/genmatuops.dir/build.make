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
include blitz/generate/CMakeFiles/genmatuops.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include blitz/generate/CMakeFiles/genmatuops.dir/compiler_depend.make

# Include the progress variables for this target.
include blitz/generate/CMakeFiles/genmatuops.dir/progress.make

# Include the compile flags for this target's objects.
include blitz/generate/CMakeFiles/genmatuops.dir/flags.make

blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o: blitz/generate/CMakeFiles/genmatuops.dir/flags.make
blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o: /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genmatuops.cpp
blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o: blitz/generate/CMakeFiles/genmatuops.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o -MF CMakeFiles/genmatuops.dir/genmatuops.cpp.o.d -o CMakeFiles/genmatuops.dir/genmatuops.cpp.o -c /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genmatuops.cpp

blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/genmatuops.dir/genmatuops.cpp.i"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genmatuops.cpp > CMakeFiles/genmatuops.dir/genmatuops.cpp.i

blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/genmatuops.dir/genmatuops.cpp.s"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate/genmatuops.cpp -o CMakeFiles/genmatuops.dir/genmatuops.cpp.s

# Object files for target genmatuops
genmatuops_OBJECTS = \
"CMakeFiles/genmatuops.dir/genmatuops.cpp.o"

# External object files for target genmatuops
genmatuops_EXTERNAL_OBJECTS =

blitz/generate/genmatuops: blitz/generate/CMakeFiles/genmatuops.dir/genmatuops.cpp.o
blitz/generate/genmatuops: blitz/generate/CMakeFiles/genmatuops.dir/build.make
blitz/generate/genmatuops: blitz/generate/CMakeFiles/genmatuops.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable genmatuops"
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/genmatuops.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
blitz/generate/CMakeFiles/genmatuops.dir/build: blitz/generate/genmatuops
.PHONY : blitz/generate/CMakeFiles/genmatuops.dir/build

blitz/generate/CMakeFiles/genmatuops.dir/clean:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate && $(CMAKE_COMMAND) -P CMakeFiles/genmatuops.dir/cmake_clean.cmake
.PHONY : blitz/generate/CMakeFiles/genmatuops.dir/clean

blitz/generate/CMakeFiles/genmatuops.dir/depend:
	cd /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2 /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-1.0.2/blitz/generate /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate /home/johannes/uni/msc/uzh/HS21/ahpc/project/blitz-build/blitz/generate/CMakeFiles/genmatuops.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : blitz/generate/CMakeFiles/genmatuops.dir/depend
