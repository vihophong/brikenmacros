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

# Utility rule file for G__fitF.

# Include the progress variables for this target.
include CMakeFiles/G__fitF.dir/progress.make

CMakeFiles/G__fitF: G__fitF.cxx
CMakeFiles/G__fitF: libfitF_rdict.pcm
CMakeFiles/G__fitF: libfitF.rootmap


G__fitF.cxx: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitFLinkDef.hh
G__fitF.cxx: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitF.hh
G__fitF.cxx: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitF.hh
G__fitF.cxx: /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitFLinkDef.hh
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating G__fitF.cxx, libfitF_rdict.pcm, libfitF.rootmap"
	/usr/bin/cmake3 -E env LD_LIBRARY_PATH=/home/phong/root_v6.14/lib:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:/home/phong/projects/geant4/geant4-qt/lib64:/home/phong/lib:/home/phong/root_v6.14/lib:/opt/rcdaq/install/lib:/usr/lib64 /home/phong/root_v6.14/bin/rootcling -v2 -f G__fitF.cxx -s /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/libfitF.so -rml libfitF.so -rmf /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/libfitF.rootmap -I/home/phong/root_v6.14/include -I/home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4 /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitF.hh /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4/include/fitFLinkDef.hh

libfitF_rdict.pcm: G__fitF.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libfitF_rdict.pcm

libfitF.rootmap: G__fitF.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libfitF.rootmap

G__fitF: CMakeFiles/G__fitF
G__fitF: G__fitF.cxx
G__fitF: libfitF_rdict.pcm
G__fitF: libfitF.rootmap
G__fitF: CMakeFiles/G__fitF.dir/build.make

.PHONY : G__fitF

# Rule to build all files generated by this target.
CMakeFiles/G__fitF.dir/build: G__fitF

.PHONY : CMakeFiles/G__fitF.dir/build

CMakeFiles/G__fitF.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/G__fitF.dir/cmake_clean.cmake
.PHONY : CMakeFiles/G__fitF.dir/clean

CMakeFiles/G__fitF.dir/depend:
	cd /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4 /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/pgenfitv4 /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit /home/phong/briken17/brikenmacrosupdate/brikenmacros/fitting/pgenfit/benchmarkPgenfitv4/build-penfit/CMakeFiles/G__fitF.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/G__fitF.dir/depend
