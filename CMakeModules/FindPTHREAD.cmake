#############################################################################
#
# $Id: FindPTHREAD.cmake 3405 2011-10-17 08:13:07Z fspindle $
#
# This file is part of the ViSP software.
# Copyright (C) 2005 - 2011 by INRIA. All rights reserved.
# 
# This software is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# ("GPL") version 2 as published by the Free Software Foundation.
# See the file LICENSE.txt at the root directory of this source
# distribution for additional information about the GNU GPL.
#
# For using ViSP with software that can not be combined with the GNU
# GPL, please contact INRIA about acquiring a ViSP Professional 
# Edition License.
#
# See http://www.irisa.fr/lagadic/visp/visp.html for more information.
# 
# This software was developed at:
# INRIA Rennes - Bretagne Atlantique
# Campus Universitaire de Beaulieu
# 35042 Rennes Cedex
# France
# http://www.irisa.fr/lagadic
#
# If you have questions regarding the use of this file, please contact
# INRIA at visp@inria.fr
# 
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# Try to find pthread library.
# Once run this will define: 
#
# PTHREAD_FOUND
# PTHREAD_INCLUDE_DIRS
# PTHREAD_LIBRARIES
#
# Authors:
# Fabien Spindler
#
#############################################################################

#IF(NOT UNIX AND NOT WIN32)
#  SET(PTHREAD_FOUND FALSE)
#ELSE(NOT UNIX AND NOT WIN32)
  
  FIND_PATH(PTHREAD_INCLUDE_DIR pthread.h
    /usr/include
    "$ENV{PTHREAD_INCLUDE_DIR}"
    "$ENV{PTHREAD_HOME}/include"
	)
  #MESSAGE("DBG PTHREAD_INCLUDE_DIR=${PTHREAD_INCLUDE_DIR}")  
  
  # pthreadVSE pthreadGCE pthreadGC pthreadVC1 pthreadVC2 are comming from web
  FIND_LIBRARY(PTHREAD_LIBRARY
    NAMES pthread pthreadGC2 pthreadVSE pthreadGCE pthreadGC pthreadVC1 pthreadVC2
    PATHS
    /usr/lib
    /usr/local/lib
    /lib    
    "$ENV{PTHREAD_LIBRARY_DIR}"
    "$ENV{PTHREAD_HOME}/lib"
    )

  #MESSAGE(STATUS "DBG PTHREAD_LIBRARY=${PTHREAD_LIBRARY}")
  
  ## --------------------------------
  
  IF(PTHREAD_LIBRARY)
    SET(PTHREAD_LIBRARIES ${PTHREAD_LIBRARY})
  ELSE(PTHREAD_LIBRARY)
    #MESSAGE(SEND_ERROR "pthread library not found.")
  ENDIF(PTHREAD_LIBRARY)
  
  IF(NOT PTHREAD_INCLUDE_DIR)
    #MESSAGE(SEND_ERROR "pthread include dir not found.")
  ENDIF(NOT PTHREAD_INCLUDE_DIR)
  
  IF(PTHREAD_LIBRARIES AND PTHREAD_INCLUDE_DIR)
    SET(PTHREAD_INCLUDE_DIRS ${PTHREAD_INCLUDE_DIR})
    SET(PTHREAD_FOUND TRUE)
  ELSE(PTHREAD_LIBRARIES AND PTHREAD_INCLUDE_DIR)
    SET(PTHREAD_FOUND FALSE)
  ENDIF(PTHREAD_LIBRARIES AND PTHREAD_INCLUDE_DIR)
  
  MARK_AS_ADVANCED(
    PTHREAD_INCLUDE_DIR
    PTHREAD_LIBRARY
  )
  #MESSAGE(STATUS "PTHREAD_FOUND : ${PTHREAD_FOUND}")

#ENDIF(NOT UNIX AND NOT WIN32)
