# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests"

# Include any dependencies generated for this target.
include src/CMakeFiles/file_writer.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/file_writer.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/file_writer.dir/flags.make

src/CMakeFiles/file_writer.dir/file_writer.cpp.o: src/CMakeFiles/file_writer.dir/flags.make
src/CMakeFiles/file_writer.dir/file_writer.cpp.o: ../src/file_writer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/file_writer.dir/file_writer.cpp.o"
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/file_writer.dir/file_writer.cpp.o -c "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/src/file_writer.cpp"

src/CMakeFiles/file_writer.dir/file_writer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/file_writer.dir/file_writer.cpp.i"
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/src/file_writer.cpp" > CMakeFiles/file_writer.dir/file_writer.cpp.i

src/CMakeFiles/file_writer.dir/file_writer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/file_writer.dir/file_writer.cpp.s"
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/src/file_writer.cpp" -o CMakeFiles/file_writer.dir/file_writer.cpp.s

src/CMakeFiles/file_writer.dir/file_writer.cpp.o.requires:

.PHONY : src/CMakeFiles/file_writer.dir/file_writer.cpp.o.requires

src/CMakeFiles/file_writer.dir/file_writer.cpp.o.provides: src/CMakeFiles/file_writer.dir/file_writer.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/file_writer.dir/build.make src/CMakeFiles/file_writer.dir/file_writer.cpp.o.provides.build
.PHONY : src/CMakeFiles/file_writer.dir/file_writer.cpp.o.provides

src/CMakeFiles/file_writer.dir/file_writer.cpp.o.provides.build: src/CMakeFiles/file_writer.dir/file_writer.cpp.o


# Object files for target file_writer
file_writer_OBJECTS = \
"CMakeFiles/file_writer.dir/file_writer.cpp.o"

# External object files for target file_writer
file_writer_EXTERNAL_OBJECTS =

src/libfile_writer.so: src/CMakeFiles/file_writer.dir/file_writer.cpp.o
src/libfile_writer.so: src/CMakeFiles/file_writer.dir/build.make
src/libfile_writer.so: src/CMakeFiles/file_writer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libfile_writer.so"
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/file_writer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/file_writer.dir/build: src/libfile_writer.so

.PHONY : src/CMakeFiles/file_writer.dir/build

src/CMakeFiles/file_writer.dir/requires: src/CMakeFiles/file_writer.dir/file_writer.cpp.o.requires

.PHONY : src/CMakeFiles/file_writer.dir/requires

src/CMakeFiles/file_writer.dir/clean:
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" && $(CMAKE_COMMAND) -P CMakeFiles/file_writer.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/file_writer.dir/clean

src/CMakeFiles/file_writer.dir/depend:
	cd "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar" "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/src" "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests" "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src" "/mnt/c/Users/gc2016/OneDrive - Imperial College London/ACSE/ACSE-4.3/acse-4-sph-morar/tests/src/CMakeFiles/file_writer.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : src/CMakeFiles/file_writer.dir/depend

