# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/brendan/Projects/Textiles

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/brendan/Projects/Textiles

# Include any dependencies generated for this target.
include CMakeFiles/ThreeByThreeA.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ThreeByThreeA.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ThreeByThreeA.dir/flags.make

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o: CMakeFiles/ThreeByThreeA.dir/flags.make
CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o: ThreeByThreeA.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/brendan/Projects/Textiles/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o -c /home/brendan/Projects/Textiles/ThreeByThreeA.cpp

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/brendan/Projects/Textiles/ThreeByThreeA.cpp > CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.i

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/brendan/Projects/Textiles/ThreeByThreeA.cpp -o CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.s

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.requires:
.PHONY : CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.requires

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.provides: CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.requires
	$(MAKE) -f CMakeFiles/ThreeByThreeA.dir/build.make CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.provides.build
.PHONY : CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.provides

CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.provides.build: CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o

# Object files for target ThreeByThreeA
ThreeByThreeA_OBJECTS = \
"CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o"

# External object files for target ThreeByThreeA
ThreeByThreeA_EXTERNAL_OBJECTS =

ThreeByThreeA: CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o
ThreeByThreeA: CMakeFiles/ThreeByThreeA.dir/build.make
ThreeByThreeA: liboutput.a
ThreeByThreeA: libtextileHelper.a
ThreeByThreeA: CMakeFiles/ThreeByThreeA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ThreeByThreeA"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ThreeByThreeA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ThreeByThreeA.dir/build: ThreeByThreeA
.PHONY : CMakeFiles/ThreeByThreeA.dir/build

CMakeFiles/ThreeByThreeA.dir/requires: CMakeFiles/ThreeByThreeA.dir/ThreeByThreeA.cpp.o.requires
.PHONY : CMakeFiles/ThreeByThreeA.dir/requires

CMakeFiles/ThreeByThreeA.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ThreeByThreeA.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ThreeByThreeA.dir/clean

CMakeFiles/ThreeByThreeA.dir/depend:
	cd /home/brendan/Projects/Textiles && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/brendan/Projects/Textiles /home/brendan/Projects/Textiles /home/brendan/Projects/Textiles /home/brendan/Projects/Textiles /home/brendan/Projects/Textiles/CMakeFiles/ThreeByThreeA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ThreeByThreeA.dir/depend

