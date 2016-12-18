# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_SOURCE_DIR = /home/weiguozheng/myapps/project_map

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/weiguozheng/myapps/project_map

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip
.PHONY : install/strip/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/weiguozheng/myapps/project_map/CMakeFiles /home/weiguozheng/myapps/project_map/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/weiguozheng/myapps/project_map/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named map

# Build rule for target.
map: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 map
.PHONY : map

# fast build rule for target.
map/fast:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/build
.PHONY : map/fast

src/constant.o: src/constant.cpp.o
.PHONY : src/constant.o

# target to build an object file
src/constant.cpp.o:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/constant.cpp.o
.PHONY : src/constant.cpp.o

src/constant.i: src/constant.cpp.i
.PHONY : src/constant.i

# target to preprocess a source file
src/constant.cpp.i:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/constant.cpp.i
.PHONY : src/constant.cpp.i

src/constant.s: src/constant.cpp.s
.PHONY : src/constant.s

# target to generate assembly for a file
src/constant.cpp.s:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/constant.cpp.s
.PHONY : src/constant.cpp.s

src/gene.o: src/gene.cpp.o
.PHONY : src/gene.o

# target to build an object file
src/gene.cpp.o:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/gene.cpp.o
.PHONY : src/gene.cpp.o

src/gene.i: src/gene.cpp.i
.PHONY : src/gene.i

# target to preprocess a source file
src/gene.cpp.i:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/gene.cpp.i
.PHONY : src/gene.cpp.i

src/gene.s: src/gene.cpp.s
.PHONY : src/gene.s

# target to generate assembly for a file
src/gene.cpp.s:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/gene.cpp.s
.PHONY : src/gene.cpp.s

src/main.o: src/main.cpp.o
.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

src/map.o: src/map.cpp.o
.PHONY : src/map.o

# target to build an object file
src/map.cpp.o:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/map.cpp.o
.PHONY : src/map.cpp.o

src/map.i: src/map.cpp.i
.PHONY : src/map.i

# target to preprocess a source file
src/map.cpp.i:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/map.cpp.i
.PHONY : src/map.cpp.i

src/map.s: src/map.cpp.s
.PHONY : src/map.s

# target to generate assembly for a file
src/map.cpp.s:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/map.cpp.s
.PHONY : src/map.cpp.s

src/mole.o: src/mole.cpp.o
.PHONY : src/mole.o

# target to build an object file
src/mole.cpp.o:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/mole.cpp.o
.PHONY : src/mole.cpp.o

src/mole.i: src/mole.cpp.i
.PHONY : src/mole.i

# target to preprocess a source file
src/mole.cpp.i:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/mole.cpp.i
.PHONY : src/mole.cpp.i

src/mole.s: src/mole.cpp.s
.PHONY : src/mole.s

# target to generate assembly for a file
src/mole.cpp.s:
	$(MAKE) -f CMakeFiles/map.dir/build.make CMakeFiles/map.dir/src/mole.cpp.s
.PHONY : src/mole.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... map"
	@echo "... rebuild_cache"
	@echo "... src/constant.o"
	@echo "... src/constant.i"
	@echo "... src/constant.s"
	@echo "... src/gene.o"
	@echo "... src/gene.i"
	@echo "... src/gene.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
	@echo "... src/map.o"
	@echo "... src/map.i"
	@echo "... src/map.s"
	@echo "... src/mole.o"
	@echo "... src/mole.i"
	@echo "... src/mole.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

