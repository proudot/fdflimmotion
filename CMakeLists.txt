# File : ---- CMakeLists.txt
# Date : 23/12/2011
# Author : Tristan Lecorgne
# INRIA Serpico team
#
# Description : This file can be used in order to prepare the build of a cimg project
#
# Usage :
#    - Read all the comments
#    - Put your preferences in the line that ends with a "# EDIT"
#

#-----------------------------------------------------------------------------
### Paths

#
# If your system use atypical developpement folders, uncomment and set the following variables.
# This is useful on Windows systems in particular.
# Note: \ needs to be escaped (example C:\\include)
#

#set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "your include path") # EDIT
#set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "your library path") # EDIT

#-----------------------------------------------------------------------------
### Check CMake version
IF(UNIX)
    CMAKE_MINIMUM_REQUIRED (VERSION 2.4.0 FATAL_ERROR)
    MESSAGE(STATUS "Compilation mode : UNIX")
ELSE(UNIX)
    IF(WIN32)
        CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6 FATAL_ERROR)
        MESSAGE(STATUS "Compilation mode : Win32")
    ELSE(WIN32)
        IF(APPLE)
            CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)
            IF(${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4 AND ${CMAKE_PATCH_VERSION} LESS 7)
                MESSAGE("Warning: A critical CMake bug exists in 2.4.6 and below. Trying to build Universal Binaries will result in a compile error that seems unrelated. Either avoid building Universal Binaries by changing the CMAKE_OSX_ARCHITECTURES field to list only your architecture, or upgrade to the current CVS version of CMake or a newer stable version if it exists.")
            ENDIF(${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4 AND ${CMAKE_PATCH_VERSION} LESS 7)
            MESSAGE(STATUS "Compilation mode : MAC OSX")
        ELSE(APPLE)
            CMAKE_MINIMUM_REQUIRED(VERSION 2.4.0 FATAL_ERROR)
        ENDIF(APPLE)
    ENDIF(WIN32)
ENDIF(UNIX)

#-----------------------------------------------------------------------------
### Project name
project(ISBICHALLENGE12) # EDIT : put your project name

#-----------------------------------------------------------------------------
set (PROJECT_BINARY_DIR "")
set (PROJECT_SOURCE_DIR "./")
### Include project directories.
include_directories (
)


##-----------------------------------------------------------------------------
### Add our own CMake modules to the search path
list(APPEND CMAKE_MODULE_PATH "${ISBICHALLENGE12_SOURCE_DIR}/CMakeModules") # EDIT : change with your project name

#-----------------------------------------------------------------------------
### Compilation.
option (COMPILATION_WARNING "Flags to print warning while compilation. (-Wall -W -ansi -pedantic)" OFF) # EDIT : set to ON to print warning
if (COMPILATION_WARNING)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -W -ansi -pedantic")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -W -ansi -pedantic")
else (COMPILATION_WARNING)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif (COMPILATION_WARNING)

option (COMPILATION_32BITS "Compile for 32 bits system (should be used on 64bits system having 32 bits libs installed)" OFF) # EDIT : set to ON to compile for 32 bits system
if (COMPILATION_32BITS)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
endif (COMPILATION_32BITS)

option (COMPILATION_OPTIMISATION "Flags to enable code optimization by the compiler." ON) # EDIT : set to OFF to disable code optimisation
if (COMPILATION_OPTIMISATION)
    if (MINGW OR CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    endif (MINGW OR CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		
    if (MSVC OR MSVC_IDE OR MSVC60 OR  MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
    endif (MSVC OR MSVC_IDE OR MSVC60 OR  MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10)

endif (COMPILATION_OPTIMISATION)

option (CIMG_DEBUG_FLAGS "Flags to enable code debugging. (-g -Dcimg_debug=3)" OFF) # EDIT
if (CIMG_DEBUG_FLAGS)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Dcimg_debug=3")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Dcimg_debug=3")
endif (CIMG_DEBUG_FLAGS)

option (PROFILING_FLAGS "Flags to enable code profiling. (-pg)" OFF) # EDIT
if (PROFILING_FLAGS)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -pg")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -pg")
	endif (PROFILING_FLAGS)
#-----------------------------------------------------------------------------
### Miscellaneous
option (CIMG_VT100_COLOR_FLAGS "Enable color output messages. (requires a VT100 compatible terminal)" ON)
if (CIMG_VT100_COLOR_FLAGS)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dcimg_use_vt100")
    set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Dcimg_use_vt100")
endif (CIMG_VT100_COLOR_FLAGS)

option (EXECUTABLE_IN_BUILD_ROOT "Put all executable in build root. OFF = executable in respective build folders" ON) # EDIT
if (EXECUTABLE_IN_BUILD_ROOT)
    set (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
endif (EXECUTABLE_IN_BUILD_ROOT)

#-----------------------------------------------------------------------------
### CIMG Support

#
# CIMG is very flexible and the use of libraries is often not mandatory
# But if your application absolutely require the use of a specific library and the
# user does not have it, you can crash the compilation and tell him to install it.
#
option (CIMG_DISPLAY_REQUIRED "If you don't handle cimg_display=0" OFF) # EDIT
option (CIMG_PNG_REQUIRED "If you don't handle cimg_use_png undefined" OFF) # EDIT
option (CIMG_JPEG_REQUIRED "If you don't handle cimg_use_jpeg undefined" OFF) # EDIT
option (CIMG_TIFF_REQUIRED "If you don't handle cimg_use_tiff undefined" OFF) # EDIT
option (CIMG_LAPACK_REQUIRED "If you don't handle cimg_use_lapack undefined" OFF) # EDIT
option (CIMG_BOARD_REQUIRED "If you don't handle cimg_use_board undefined" OFF) # EDIT
option (CIMG_FFMPEG_REQUIRED "If you don't handle cimg_use_ffmpeg undefined" OFF) # EDIT
option (CIMG_FFTW3_REQUIRED "If you don't handle cimg_use_fftw3 undefined" OFF) # EDIT
option (CIMG_MAGICK_REQUIRED "If you don't handle cimg_use_magick undefined" OFF) # EDIT


#
# If you don't want your application to be linked against all cimg optional libraries,
# you can disable some of them here. It is useful for releases.
#
option (USE_DISPLAY_LIBS "Set to off if your application does not use any display (=> cimg_diplay=0)" ON) # EDIT
option (USE_JPEG_LIBS "Set to off if you don't need libjpeg (=> cimg_use_jpeg undefined)" ON) # EDIT
option (USE_PNG_LIBS "Set to off if you don't need libpng (=> cimg_use_png undefined)" ON) # EDIT
option (USE_TIFF_LIBS "Set to off if you don't need libtiff (=> cimg_use_tiff undefined)" ON) # EDIT
option (USE_LAPACK_LIBS "Set to off if you don't need lapack libraries (=> cimg_use_lapack undefined)" OFF) # EDIT
option (USE_OPENMP_OPT "Set to off if you don't need openmp optimisation (=> cimg_use_openmp undefined)" OFF) # EDIT
option (USE_BOARD_LIBS "Set to off if you don't need libboard (=> cimg_use_board undefined)" OFF) # EDIT
option (USE_FFMPEG_LIBS "Set to off if you don't need libffmpeg (=> cimg_use_ffmpeg undefined)" OFF) # EDIT
option (USE_FFTW3_LIBS "Set to off if you don't need libfftw3 (=> cimg_use_fftw3 undefined)" OFF) # EDIT
option (USE_MAGICK_LIBS "Set to off if you don't need libmagick++ (=> cimg_use_magick undefined)" OFF) # EDIT

#
# Call findCIMG.cmake and use the variables defined by it
#
find_package (CIMG)
add_definitions (${CIMG_DEFINITIONS})
include_directories (${CIMG_INCLUDE_DIRS})
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CIMG_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CIMG_CXX_FLAGS}")



# ISBICHALLENGE12 specific paths
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
ADD_DEFINITIONS(-Dcimg_plugin1=\"\\\"noise_analysis.h\\\"\")

#-----------------------------------------------------------------------------
### Executables

#
# Here you can declare what executable you want to compile with cmake commands.
# For example :
# - add_executable (myExecutable myExecutable.cpp myExecutable.h otherRequiredFile.h ...)
# - add_library (myLibrary myLibrary.cpp ...)
# - target_link_libraries (myExecutable myLibrary  ...)
#
# Then, link ${CIMG_LIBRARIES} with your executables :
#
# - target_link_libraries (myExecutable ${CIMG_LIBRARIES})
#
# For more information on these commands and on other possible commands, go to cmake.org
#
# -- NOTE --
# It is strongly recommanded not to make recursive CMakeLists.txt.
# So if you need to compile myExecutable with files located in ${MyCimgProject_SOURCE_DIR}/src write :
# add_executable (myExecutable src/myExecutable.cpp src/myExecutable.h src/otherRequiredFile.h ...)
#

# Maintained binaries
ADD_EXECUTABLE(fdflimmotion ${PROJECT_SOURCE_DIR}fdflimmotion.cpp  ${PROJECT_SOURCE_DIR}sin_fit.cpp ${PROJECT_SOURCE_DIR}tinyxml2.cpp)
TARGET_LINK_LIBRARIES(fdflimmotion ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/fdflimmotion DESTINATION bin)

ADD_EXECUTABLE(FDFLIMRepair ${PROJECT_SOURCE_DIR}FDFLIMRepair.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(FDFLIMRepair ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/FDFLIMRepair DESTINATION bin)

ADD_EXECUTABLE(fourier_lifetime ${PROJECT_SOURCE_DIR}fourier_lifetime.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(fourier_lifetime ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/fourier_lifetime DESTINATION bin)

# Figures and test script
ADD_EXECUTABLE(abber-impact-tau ${PROJECT_SOURCE_DIR}abber-impact-tau.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(abber-impact-tau ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/abber-impact-tau DESTINATION bin)

ADD_EXECUTABLE(phase-iccd-photon-count ${PROJECT_SOURCE_DIR}phase-iccd-photon-count.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(phase-iccd-photon-count ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/phase-iccd-photon-count DESTINATION bin)

ADD_EXECUTABLE(phase-abber-iccd-photon-count-moving ${PROJECT_SOURCE_DIR}phase-abber-iccd-photon-count-moving.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(phase-abber-iccd-photon-count-moving ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/phase-abber-iccd-photon-count-moving DESTINATION bin)

ADD_EXECUTABLE(energy-minimum ${PROJECT_SOURCE_DIR}energy-minimum.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(energy-minimum ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/energy-minimum/ DESTINATION bin)


ADD_EXECUTABLE(lifetime ${PROJECT_SOURCE_DIR}lifetime.cpp ${PROJECT_SOURCE_DIR}sin_fit.cpp)
TARGET_LINK_LIBRARIES(lifetime ${CIMG_LIBRARIES})
INSTALL(PROGRAMS build/lifetime DESTINATION bin)
