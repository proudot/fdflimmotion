# File : ---- FindCIMG.cmake
# Date : 23/12/2011
# Author : Tristan Lecorgne
# INRIA Serpico team
#
# Description : Try to find the cimg optional libs on the system and set the following variables :
#
#        CIMG_LIBRARIES       The .so/.a/.dll/... to link against
#        CIMG_INCLUDE_DIRS    The headers directory of the CIMG_LIBRARIES
#        CIMG_DEFINITIONS     The -Dcimg_use_xxxx options
#        CIMG_C_FLAGS         The flags required to compile c files with cimg
#        CIMG_CXX_FLAGS       The flags required to compile c++ files with cimg
#
# Usage : 
#
#    - add the following lines to your main CMakeLists.txt:
#    -- set the CIMG_XXX_REQUIRED options to specify which libraries you absolutely need:
#    ---- option (CIMG_DISPLAY_REQUIRED "If you don't handle cimg_display=0" OFF)
#    ---- option (CIMG_PNG_REQUIRED "If you don't handle cimg_use_png undefined" OFF)
#    ---- option (CIMG_JPEG_REQUIRED "If you don't handle cimg_use_jpeg undefined" OFF)
#    ---- option (CIMG_TIFF_REQUIRED "If you don't handle cimg_use_tiff undefined" OFF)
#    ---- option (CIMG_LAPACK_REQUIRED "If you don't handle cimg_use_lapack undefined" OFF)
#    ---- option (CIMG_BOARD_REQUIRED "If you don't handle cimg_use_board undefined" OFF)
#    ---- option (CIMG_FFMPEG_REQUIRED "If you don't handle cimg_use_ffmpeg undefined" OFF)
#    ---- option (CIMG_FFTW3_REQUIRED "If you don't handle cimg_use_fftw3 undefined" OFF)
#    ---- option (CIMG_MAGICK_REQUIRED "If you don't handle cimg_use_magick undefined" OFF)
#    -- set the CIMG_USE_XXX options to specify which libraries you don't want to use:
#    ---- option (USE_DISPLAY_LIBS "Set to off if your application does not use any display (=> cimg_diplay=0)" ON) 
#    ---- option (USE_JPEG_LIBS "Set to off if you don't need libjpeg (=> cimg_use_jpeg undefined)" ON)
#    ---- option (USE_PNG_LIBS "Set to off if you don't need libpng (=> cimg_use_png undefined)" ON)
#    ---- option (USE_TIFF_LIBS "Set to off if you don't need libtiff (=> cimg_use_tiff undefined)" ON)
#    ---- option (USE_LAPACK_LIBS "Set to off if you don't need lapack libraries (=> cimg_use_lapack undefined)" ON)
#    ---- option (USE_OPENMP_OPT "Set to off if you don't need openmp optimisation (=> cimg_use_openmp undefined)" ON)
#    ---- option (USE_BOARD_LIBS "Set to off if you don't need libboard (=> cimg_use_board undefined)" ON)
#    ---- option (USE_FFMPEG_LIBS "Set to off if you don't need libffmpeg (=> cimg_use_ffmpeg undefined)" ON)
#    ---- option (USE_FFTW3_LIBS "Set to off if you don't need libfftw3 (=> cimg_use_fftw3 undefined)" ON)
#    ---- option (USE_MAGICK_LIBS "Set to off if you don't need libmagick++ (=> cimg_use_magick undefined)" ON)
#    -- add the CMakeModules directory to the module path:
#    ---- list(APPEND CMAKE_MODULE_PATH ${<projectName>_SOURCE_DIR}/CMakeModules)
#    -- call the execution of this file
#    ---- find_PACKAGE(CIMG)
#    -- use the variables defined by this file:
#    ---- add_definitions (${CIMG_DEFINITIONS})
#    ---- include_directories (${CIMG_INCLUDE_DIRS})
#    ---- set (${CMAKE_C_FLAGS} "${CMAKE_C_FLAGS} ${CIMG_C_FLAGS}")
#    ---- set (${CMAKE_CXX_FLAGS} "${CMAKE_CXX_FLAGS} ${CIMG_CXX_FLAGS}")
#    -- for each executable that depend on Cimg.h, add the line:
#    ---- target_link_libraries (my_executable ${CIMG_LIBRARIES})
#
# Details :
#    - Find the appropriate display library according to the system (x11 on unix, gdi on windows).
#        If the library is X11, it search for pthread, which is required by cimg when using x11
#        and the xshm and xrandr optional libraries. 
#
#    - Try to find libpng, libjpeg and libtiff. Those libraries are optional but if ommited, cimg
#        will use imagemagick programs to perform load and save operations. Imagemagick is often 
#        less efficient and can't handle some specific format (for example 3D tiff images).
#
#    - Try to see if OpenMP support is possible and add the c and c++ flags.
#    
#    - Try to find the lapack, board, fftw3, magick++ and ffmpeg.
# 
#    - Finaly, check that all the libraries required by the developper are found.
#


### DISPLAY :: X11 on unix-based system and GDI on windows ###
if (USE_DISPLAY_LIBS)
    if (UNIX OR APPLE)
        find_package (X11 QUIET) # xshm xrandr are detected as well
        if (X11_FOUND)
            set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${X11_INCLUDE_DIR} )
            set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${X11_LIBRARIES})
            message (STATUS "FindCIMG.cmake : X11 found.")
            
    ### X11 extension :: XSHM ###
            if (X11_XShm_FOUND)
                set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_xshm)
                message(STATUS "FindCIMG.cmake : xshm found")
            else (X11_XShm_FOUND)
                message(STATUS "!!! FindCIMG.cmake !!! xshm NOT found. This is an optional X11 extension that improve memory management. It is recommended. If you are the developper, be sure that you handle the case where cimg_use_xshm is not defined.")
            endif (X11_XShm_FOUND)
    
    ### X11 extension :: XRANDR ###
            if (X11_Xrandr_FOUND)
                set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_xrandr)
                set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${X11_Xrandr_LIB})
                message(STATUS "FindCIMG.cmake : xrandr found")
            else (X11_Xrandr_FOUND)
                message(STATUS "!!! FindCIMG.cmake !!! xrandr NOT found")
            endif (X11_Xrandr_FOUND)
            
    ### PThread is required when using X11 display engine ###
            find_package (PTHREAD QUIET)
            if (PTHREAD_FOUND)
                set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${PTHREAD_INCLUDE_DIRS})
                set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${PTHREAD_LIBRARIES})
                message(STATUS "FindCIMG.cmake : pthread found")
            else (PTHREAD_FOUND)
                message(STATUS "!!! FindCIMG.cmake !!! pthread NOT found. pthread required by cimg for running X11. The application will have NO DISPLAY. If you are the developper, make sure that you handle the case where cimg_display is set to 0")
                set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_display=0)  
            endif (PTHREAD_FOUND)
            
        else (X11_FOUND)
            message (STATUS "!!! FindCIMG.cmake !!! X11 NOT found. The application will have NO DISPLAY. If you are the developper, make sure that you handle the case where cimg_display is set to 0")
            set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_display=0)  
            set (CIMG_NO_DISPLAY)
        endif (X11_FOUND)
    else (UNIX OR APPLE)
        if (WIN32)
            find_package(GDI QUIET)
            if (GDI_FOUND)
                set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${GDI_LIBRARIES})
                message(STATUS "FindCIMG.cmake : GDI found")
            else (GDI_FOUND)
                message(STATUS "!!! FindCIMG.cmake !!! GDI NOT found. The application will have NO DISPLAY. If you are the developper, make sure that you handle the case where cimg_display is set to 0")
                set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_display=0)  
                set (CIMG_NO_DISPLAY)
            endif (GDI_FOUND)
        endif (WIN32)
    endif (UNIX OR APPLE)
else (USE_DISPLAY_LIBS)
    message(STATUS "FindCIMG.cmake : display disabled")
    set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_display=0)    
endif (USE_DISPLAY_LIBS)

### BETTER IMAGE SUPPORT :: PNG ###
if (USE_PNG_LIBS)
    find_package (PNG QUIET) # zlib inside
    if (PNG_FOUND)
        # if png is found, then zlib is included in png_include_dir and png_libraries
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${PNG_INCLUDE_DIR})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${PNG_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_png -Dcimg_use_zlib ${PNG_DEFINITIONS})    
        message(STATUS "FindCIMG.cmake : zlib found")
        message(STATUS "FindCIMG.cmake : png found")
    else (PNG_FOUND)
        if (ZLIB_FOUND)
            set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR})
            set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${ZLIB_LIBRARIES})
            set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_zlib) 
            message(STATUS "FindCIMG.cmake : zlib found")
        else (ZLIB_FOUND)
            message(STATUS "!!! FindCIMG.cmake !!! zlib NOT found. zlib is required by libpng.")
        endif(ZLIB_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! png NOT found. You may not be able to perform some png image manipulation... If you are the developper, be sure you handle the case where cimg_use_png is not defined")
    endif (PNG_FOUND)
endif (USE_PNG_LIBS)
   
### BETTER IMAGE SUPPORT :: JPEG ###
if (USE_JPEG_LIBS)
    find_package (JPEG QUIET)
    if (JPEG_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${JPEG_INCLUDE_DIR})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${JPEG_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_jpeg)    
        message(STATUS "FindCIMG.cmake : jpeg found")
    else (JPEG_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! jpeg NOT found. You may not be able to perform some jpeg image manipulation... If you are the developper, be sure you handle the case where cimg_use_jpeg is not defined")
    endif (JPEG_FOUND)
endif (USE_JPEG_LIBS)

### BETTER IMAGE SUPPORT :: TIFF ###
if (USE_TIFF_LIBS)
    find_package (TIFF QUIET)
    if (TIFF_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${TIFF_INCLUDE_DIR})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${TIFF_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_tiff)    
        message(STATUS "FindCIMG.cmake : tiff found")
    else (TIFF_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! tiff NOT found. You may not be able to perform some tiff image manipulation... If you are the developper, be sure you handle the case where cimg_use_tiff is not defined")
    endif (TIFF_FOUND)
endif (USE_TIFF_LIBS)

### OPTIMISATION :: OPENMP enable multiprocessor computing support ###
if (USE_OPENMP_OPT)
    find_package (OpenMP QUIET)
    if (OPENMP_FOUND)
        set (CIMG_C_FLAGS "${CIMG_C_FLAGS}  ${OpenMP_C_FLAGS}")
        set (CIMG_CXX_FLAGS "${CIMG_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_openmp)    
        message(STATUS "FindCIMG.cmake : openmp found")
    else (OPENMP_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! openmp NOT found. The program will not use multi-processor optimisations ... If you are the developper, be sure you handle the case where cimg_use_openMP is not defined")
    endif (OPENMP_FOUND)
endif (USE_OPENMP_OPT)

### LINEAR ALGEBRA PACKAGE :: LAPACK ###
if (USE_LAPACK_LIBS)
    find_package(LAPACK QUIET)
    if (LAPACK_FOUND)
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${LAPACK_LIBRARIES})
        set (CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        set (CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        set (CMAKE_MODULE_LINKER_FLAGS ${CMAKE_MODULE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_lapack)
    else (LAPACK_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! lapack NOT found. You may not be able to perform some linear algebric computation... If you are the developper, be sure you handle the case where cimg_use_lapack is not defined")
    endif (LAPACK_FOUND)
endif (USE_LAPACK_LIBS) 


### 3D GRAPHICS :: BOARD ###
if (USE_BOARD_LIBS)
    find_package (Board QUIET)
    if (BOARD_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${BOARD_INCLUDE_DIRS})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${BOARD_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_board ${BOARD_DEFINITIONS})    
        message(STATUS "FindCIMG.cmake : board found")
    else (BOARD_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! board NOT found. You may not be able to perform some 3D graphics computation... If you are the developper, be sure you handle the case where cimg_use_board is not defined")
    endif (BOARD_FOUND)
endif (USE_BOARD_LIBS)

### VIDEO :: FFMPEG ###
if (USE_FFMPEG_LIBS)
    find_package (FFMPEG QUIET) # include avformat, avcodec
    if (FFMPEG_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${FFMPEG_INCLUDE_DIRS})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${FFMPEG_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_ffmpeg)    
        message(STATUS "FindCIMG.cmake : ffmpeg found")
    else (FFMPEG_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! ffmpeg NOT found. You may not be able to perform some video manipulation... If you are the developper, be sure you handle the case where cimg_use_ffmpeg is not defined")
    endif (FFMPEG_FOUND)
endif (USE_FFMPEG_LIBS)

### FAST FOURIER TRANSFORM :: FFTW3 ###
if (USE_FFTW3_LIBS)
    find_package (Fftw3 QUIET)
    if (FFTW3_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIR})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${FFTW3_LIBRARIES})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_fftw3)    
        message(STATUS "FindCIMG.cmake : fftw3 found")
    else (FFTW3_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! fftw3 NOT found. You may not be able to perform some fft computation... If you are the developper, be sure you handle the case where cimg_use_fftw3 is not defined")
    endif (FFTW3_FOUND)
endif (USE_FFTW3_LIBS)

### IMAGE FILE MANAGEMENT :: MAGICK++ ###
if (USE_MAGICK_LIBS)
    find_package (Magick QUIET)
    if (MAGICK_FOUND)
        set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIR})
        set (CIMG_LIBRARIES ${CIMG_LIBRARIES} ${MAGICK_LIBRARY} ${MAGICK++_LIBRARY})
        set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} -Dcimg_use_magick)    
        message(STATUS "FindCIMG.cmake : magick++ found")
    else (MAGICK_FOUND)
        message(STATUS "!!! FindCIMG.cmake !!! magick++ NOT found. You may not be able to perform some image load/save operations... If you are the developper, be sure you handle the case where cimg_use_magick is not defined")
    endif (MAGICK_FOUND)
endif (USE_MAGICK_LIBS)


### Check Required libraries ###
if (CIMG_DISPLAY_REQUIRED AND CIMG_NO_DISPLAY)
    message (FATAL_ERROR "You need a display engine such as X11 or GDI to compile this program. Please install libs and developpement headers")
endif (CIMG_DISPLAY_REQUIRED AND CIMG_NO_DISPLAY)

if (CIMG_PNG_REQUIRED AND NOT PNG_FOUND)
    message (FATAL_ERROR "You need libpng to compile this program. Please install libs and developpement headers")
endif (CIMG_PNG_REQUIRED AND NOT PNG_FOUND)

if (CIMG_JPEG_REQUIRED AND NOT JPEG_FOUND)
    message (FATAL_ERROR "You need libjpeg to compile this program. Please install libs and developpement headers")
endif (CIMG_JPEG_REQUIRED AND NOT JPEG_FOUND)

if (CIMG_TIFF_REQUIRED AND NOT TIFF_FOUND)
    message (FATAL_ERROR "You need libtiff to compile this program. Please install libs and developpement headers")
endif (CIMG_TIFF_REQUIRED AND NOT TIFF_FOUND)

if (CIMG_LAPACK_REQUIRED AND NOT LAPACK_FOUND)
    message (FATAL_ERROR "You need liblapack to compile this program. Please install libs and developpement headers")
endif (CIMG_LAPACK_REQUIRED AND NOT LAPACK_FOUND)

if (CIMG_BOARD_REQUIRED AND NOT BOARD_FOUND)
    message (FATAL_ERROR "You need libboard to compile this program. Please install libs and developpement headers")
endif (CIMG_BOARD_REQUIRED AND NOT BOARD_FOUND)

if (CIMG_FFMPEG_REQUIRED AND NOT FFMPEG_FOUND)
    message (FATAL_ERROR "You need libffmpeg to compile this program. Please install libs and developpement headers")
endif (CIMG_FFMPEG_REQUIRED AND NOT FFMPEG_FOUND)

if (CIMG_FFTW3_REQUIRED AND NOT FFTW3_FOUND)
    message (FATAL_ERROR "You need libfftw3 to compile this program. Please install libs and developpement headers")
endif (CIMG_FFTW3_REQUIRED AND NOT FFTW3_FOUND)

if (CIMG_MAGICK_REQUIRED AND NOT MAGICK_FOUND)
    message (FATAL_ERROR "You need libmagic++ to compile this program. Please install libs and developpement headers")
endif (CIMG_MAGICK_REQUIRED AND NOT MAGICK_FOUND)

set (CIMG_INCLUDE_DIRS ${CIMG_INCLUDE_DIRS} CACHE STRING "include directories for cimg dependancies")
set (CIMG_LIBRARIES ${CIMG_LIBRARIES} CACHE STRING "cimg required and optional 3rd party libraries")
set (CIMG_DEFINITIONS ${CIMG_DEFINITIONS} CACHE STRING "cimg_use_xxx defines")
set (CIMG_C_FLAGS ${CIMG_C_FLAGS}  CACHE STRING "c flags for cimg")
set (CIMG_CXX_FLAGS ${CIMG_CXX_FLAGS} CACHE STRING "c++ flags for cimg")



