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
CMAKE_COMMAND = /usr/local/Cellar/cmake/2.8.10.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/2.8.10.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/Cellar/cmake/2.8.10.2/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/brberg/Projects/Textiles

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/brberg/Projects/Textiles/build

# Include any dependencies generated for this target.
include CMakeFiles/output.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/output.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/output.dir/flags.make

CMakeFiles/output.dir/output.cpp.o: CMakeFiles/output.dir/flags.make
CMakeFiles/output.dir/output.cpp.o: ../output.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/brberg/Projects/Textiles/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/output.dir/output.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/output.dir/output.cpp.o -c /Users/brberg/Projects/Textiles/output.cpp

CMakeFiles/output.dir/output.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/output.dir/output.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/brberg/Projects/Textiles/output.cpp > CMakeFiles/output.dir/output.cpp.i

CMakeFiles/output.dir/output.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/output.dir/output.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/brberg/Projects/Textiles/output.cpp -o CMakeFiles/output.dir/output.cpp.s

CMakeFiles/output.dir/output.cpp.o.requires:
.PHONY : CMakeFiles/output.dir/output.cpp.o.requires

CMakeFiles/output.dir/output.cpp.o.provides: CMakeFiles/output.dir/output.cpp.o.requires
	$(MAKE) -f CMakeFiles/output.dir/build.make CMakeFiles/output.dir/output.cpp.o.provides.build
.PHONY : CMakeFiles/output.dir/output.cpp.o.provides

CMakeFiles/output.dir/output.cpp.o.provides.build: CMakeFiles/output.dir/output.cpp.o

# Object files for target output
output_OBJECTS = \
"CMakeFiles/output.dir/output.cpp.o"

# External object files for target output
output_EXTERNAL_OBJECTS =

liboutput.a: CMakeFiles/output.dir/output.cpp.o
liboutput.a: CMakeFiles/output.dir/build.make
liboutput.a: CMakeFiles/output.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library liboutput.a"
	$(CMAKE_COMMAND) -P CMakeFiles/output.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/output.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/output.dir/build: liboutput.a
.PHONY : CMakeFiles/output.dir/build

CMakeFiles/output.dir/requires: CMakeFiles/output.dir/output.cpp.o.requires
.PHONY : CMakeFiles/output.dir/requires

CMakeFiles/output.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/output.dir/cmake_clean.cmake
.PHONY : CMakeFiles/output.dir/clean

CMakeFiles/output.dir/depend:
	cd /Users/brberg/Projects/Textiles/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/brberg/Projects/Textiles /Users/brberg/Projects/Textiles /Users/brberg/Projects/Textiles/build /Users/brberg/Projects/Textiles/build /Users/brberg/Projects/Textiles/build/CMakeFiles/output.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/output.dir/depend

