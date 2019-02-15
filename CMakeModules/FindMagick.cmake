# - Find Magick
# This module finds if the Magick++ libraries and headers (distributed
# with ImageMagick are installed.  If you are looking for the ImageMagick
# programs please see FIND_PACKAGE(ImageMagick).
#
# This module sets the following variables:
#
# MAGICK_FOUND
#    True if the Magick++ library was found
# MAGICK_INCLUDE_DIR
#    The include path of the Magick++.h header file
# MAGICK_LIBRARY
#    The location of the Magick library
# MAGICK++_LIBRARY
#    THe location of the Magick++ library
#

# Find the Magick++ include path
FIND_PATH(MAGICK_INCLUDE_DIR Magick++.h
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ImageMagick\\Current;BinPath]/include"
)

FIND_LIBRARY(MAGICK_LIBRARY Magick CORE_RL_Magick_
    PATHS 
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ImageMagick\\Current;BinPath]/lib"
)

FIND_LIBRARY(MAGICK++_LIBRARY Magick++ CORE_RL_Magick++_
    PATHS
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ImageMagick\\Current;BinPath]/lib"
)

IF(MAGICK_INCLUDE_DIR AND MAGICK_LIBRARY AND MAGICK++_LIBRARY)
    SET(MAGICK_FOUND true)
    SET(MAGICK_INCLUDE_DIRS ${MAGICK_INCLUDE_DIR})
    SET(MAGICK_LIBRARIES    ${MAGICK_LIBRARY} ${MAGICK++_LIBRARY})
ENDIF(MAGICK_INCLUDE_DIR AND MAGICK_LIBRARY AND MAGICK++_LIBRARY)

IF(MAGICK_FOUND)
    IF(NOT MAGICK_FIND_QUIETLY)
        MESSAGE(STATUS "Found Magick: ${MAGICK_LIBRARY}")
        MESSAGE(STATUS "Found Magick++: ${MAGICK++_LIBRARY}")
    ENDIF(NOT MAGICK_FIND_QUIETLY)
ELSE(MAGICK_FOUND) 
    IF(MAGICK_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find the Magick library")
    ENDIF(MAGICK_FIND_REQUIRED)
ENDIF(MAGICK_FOUND)
