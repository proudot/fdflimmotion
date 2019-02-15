# - Try to find Board
# Once done this will define
#
#  BOARD_FOUND - system has Board
#  BOARD_INCLUDE_DIRS - the Board include directory
#  BOARD_LIBRARIES - Link these to use Board
#  BOARD_DEFINITIONS - Compiler switches required for using Board
#
#  Copyright (c) 2010 Roman Putanowicz <putanowr@l5.pk.edu.pl>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if (BOARD_LIBRARIES AND BOARD_INCLUDE_DIRS)
  # in cache already
  set(BOARD_FOUND TRUE)
else (BOARD_LIBRARIES AND BOARD_INCLUDE_DIRS)
  find_path(BOARD_INCLUDE_DIR
    NAMES
      Board.h
    PATHS
      "${BASE_DIR}/include"
      /usr/include
      /usr/local/include
      /opt/local/include
      /sw/include
  )

  find_library(BOARD_LIBRARY
    NAMES
      board
    PATHS
      "${BASE_DIR_LIB}"
      /usr/lib
      /usr/local/lib
      /opt/local/lib
      /sw/lib
  )

  set(BOARD_INCLUDE_DIRS
    ${BOARD_INCLUDE_DIR} CACHE PATH "Path to Board headers"
  )

  set(BOARD_LIBRARIES
      ${BOARD_LIBRARY} CACHE STRING "Directories to be linked to use Board"
  )

  if (BOARD_INCLUDE_DIRS AND BOARD_LIBRARIES)
     set(BOARD_FOUND TRUE)
  endif (BOARD_INCLUDE_DIRS AND BOARD_LIBRARIES)

  if (BOARD_FOUND)
    if (NOT Board_FIND_QUIETLY)
      message(STATUS "Found Board: ${BOARD_LIBRARIES}")
    endif (NOT Board_FIND_QUIETLY)
  else (BOARD_FOUND)
    if (Board_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find Board")
    endif (Board_FIND_REQUIRED)
  endif (BOARD_FOUND)

  # show the BOARD_INCLUDE_DIRS and BOARD_LIBRARIES variables only in the advanced view
  mark_as_advanced(BOARD_INCLUDE_DIRS BOARD_LIBRARIES)

endif (BOARD_LIBRARIES AND BOARD_INCLUDE_DIRS)

