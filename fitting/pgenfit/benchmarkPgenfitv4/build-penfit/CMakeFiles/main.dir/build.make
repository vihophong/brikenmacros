# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/src/decaypath.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/decaypath.cc.o: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/decaypath.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/src/decaypath.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/src/decaypath.cc.o -c /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/decaypath.cc

CMakeFiles/main.dir/src/decaypath.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/decaypath.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/decaypath.cc > CMakeFiles/main.dir/src/decaypath.cc.i

CMakeFiles/main.dir/src/decaypath.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/decaypath.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/decaypath.cc -o CMakeFiles/main.dir/src/decaypath.cc.s

CMakeFiles/main.dir/src/decaypath.cc.o.requires:

.PHONY : CMakeFiles/main.dir/src/decaypath.cc.o.requires

CMakeFiles/main.dir/src/decaypath.cc.o.provides: CMakeFiles/main.dir/src/decaypath.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/decaypath.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/src/decaypath.cc.o.provides

CMakeFiles/main.dir/src/decaypath.cc.o.provides.build: CMakeFiles/main.dir/src/decaypath.cc.o


CMakeFiles/main.dir/src/fitF.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/fitF.cc.o: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/src/fitF.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/src/fitF.cc.o -c /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF.cc

CMakeFiles/main.dir/src/fitF.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/fitF.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF.cc > CMakeFiles/main.dir/src/fitF.cc.i

CMakeFiles/main.dir/src/fitF.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/fitF.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF.cc -o CMakeFiles/main.dir/src/fitF.cc.s

CMakeFiles/main.dir/src/fitF.cc.o.requires:

.PHONY : CMakeFiles/main.dir/src/fitF.cc.o.requires

CMakeFiles/main.dir/src/fitF.cc.o.provides: CMakeFiles/main.dir/src/fitF.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/fitF.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/src/fitF.cc.o.provides

CMakeFiles/main.dir/src/fitF.cc.o.provides.build: CMakeFiles/main.dir/src/fitF.cc.o


CMakeFiles/main.dir/src/fitF_auxiliary.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/fitF_auxiliary.cc.o: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF_auxiliary.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/src/fitF_auxiliary.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/src/fitF_auxiliary.cc.o -c /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF_auxiliary.cc

CMakeFiles/main.dir/src/fitF_auxiliary.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/fitF_auxiliary.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF_auxiliary.cc > CMakeFiles/main.dir/src/fitF_auxiliary.cc.i

CMakeFiles/main.dir/src/fitF_auxiliary.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/fitF_auxiliary.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/fitF_auxiliary.cc -o CMakeFiles/main.dir/src/fitF_auxiliary.cc.s

CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.requires:

.PHONY : CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.requires

CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.provides: CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.provides

CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.provides.build: CMakeFiles/main.dir/src/fitF_auxiliary.cc.o


CMakeFiles/main.dir/src/unbinfit.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/unbinfit.cc.o: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/unbinfit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/main.dir/src/unbinfit.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/src/unbinfit.cc.o -c /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/unbinfit.cc

CMakeFiles/main.dir/src/unbinfit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/unbinfit.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/unbinfit.cc > CMakeFiles/main.dir/src/unbinfit.cc.i

CMakeFiles/main.dir/src/unbinfit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/unbinfit.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/src/unbinfit.cc -o CMakeFiles/main.dir/src/unbinfit.cc.s

CMakeFiles/main.dir/src/unbinfit.cc.o.requires:

.PHONY : CMakeFiles/main.dir/src/unbinfit.cc.o.requires

CMakeFiles/main.dir/src/unbinfit.cc.o.provides: CMakeFiles/main.dir/src/unbinfit.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/src/unbinfit.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/src/unbinfit.cc.o.provides

CMakeFiles/main.dir/src/unbinfit.cc.o.provides.build: CMakeFiles/main.dir/src/unbinfit.cc.o


CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cc.o -c /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/main.cc

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/main.cc > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/main.cc -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/main.cc.o.requires:

.PHONY : CMakeFiles/main.dir/main.cc.o.requires

CMakeFiles/main.dir/main.cc.o.provides: CMakeFiles/main.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cc.o.provides

CMakeFiles/main.dir/main.cc.o.provides.build: CMakeFiles/main.dir/main.cc.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/src/decaypath.cc.o" \
"CMakeFiles/main.dir/src/fitF.cc.o" \
"CMakeFiles/main.dir/src/fitF_auxiliary.cc.o" \
"CMakeFiles/main.dir/src/unbinfit.cc.o" \
"CMakeFiles/main.dir/main.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/src/decaypath.cc.o
main: CMakeFiles/main.dir/src/fitF.cc.o
main: CMakeFiles/main.dir/src/fitF_auxiliary.cc.o
main: CMakeFiles/main.dir/src/unbinfit.cc.o
main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/build.make
main: /home/phong/root_v6.14/lib/libCore.so
main: /home/phong/root_v6.14/lib/libImt.so
main: /home/phong/root_v6.14/lib/libRIO.so
main: /home/phong/root_v6.14/lib/libNet.so
main: /home/phong/root_v6.14/lib/libHist.so
main: /home/phong/root_v6.14/lib/libGraf.so
main: /home/phong/root_v6.14/lib/libGraf3d.so
main: /home/phong/root_v6.14/lib/libGpad.so
main: /home/phong/root_v6.14/lib/libTree.so
main: /home/phong/root_v6.14/lib/libTreePlayer.so
main: /home/phong/root_v6.14/lib/libRint.so
main: /home/phong/root_v6.14/lib/libPostscript.so
main: /home/phong/root_v6.14/lib/libMatrix.so
main: /home/phong/root_v6.14/lib/libPhysics.so
main: /home/phong/root_v6.14/lib/libMathCore.so
main: /home/phong/root_v6.14/lib/libThread.so
main: /home/phong/root_v6.14/lib/libMultiProc.so
main: /home/phong/root_v6.14/lib/libRooFit.so
main: /home/phong/root_v6.14/lib/libRooFitCore.so
main: /home/phong/root_v6.14/lib/libRooStats.so
main: libfitF.so
main: /home/phong/root_v6.14/lib/libCore.so
main: /home/phong/root_v6.14/lib/libImt.so
main: /home/phong/root_v6.14/lib/libRIO.so
main: /home/phong/root_v6.14/lib/libNet.so
main: /home/phong/root_v6.14/lib/libHist.so
main: /home/phong/root_v6.14/lib/libGraf.so
main: /home/phong/root_v6.14/lib/libGraf3d.so
main: /home/phong/root_v6.14/lib/libGpad.so
main: /home/phong/root_v6.14/lib/libTree.so
main: /home/phong/root_v6.14/lib/libTreePlayer.so
main: /home/phong/root_v6.14/lib/libRint.so
main: /home/phong/root_v6.14/lib/libPostscript.so
main: /home/phong/root_v6.14/lib/libMatrix.so
main: /home/phong/root_v6.14/lib/libPhysics.so
main: /home/phong/root_v6.14/lib/libMathCore.so
main: /home/phong/root_v6.14/lib/libThread.so
main: /home/phong/root_v6.14/lib/libMultiProc.so
main: /home/phong/root_v6.14/lib/libRooFit.so
main: /home/phong/root_v6.14/lib/libRooFitCore.so
main: /home/phong/root_v6.14/lib/libRooStats.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/src/decaypath.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/src/fitF.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/src/fitF_auxiliary.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/src/unbinfit.cc.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cc.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4 /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4 /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend
