#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC and other Celeritas Developers.
# See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
# CONFIGURE FILE
#----------------------------------------------------------------------------#

set(_CONFIG_NAME "celeritas_config.h")
configure_file("${_CONFIG_NAME}.in" "${_CONFIG_NAME}" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_CONFIG_NAME}"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

#----------------------------------------------------------------------------#
# LIBRARY
#----------------------------------------------------------------------------#
set(SOURCES
  base/KernelParamCalculator.cu
)
set(LIBRARY_DEPS
)

if(CELERITAS_USE_CUDA)
  list(APPEND SOURCES
    base/KernelParamCalculator.cu
  )
endif()

if(CELERITAS_USE_MPI)
  list(APPEND SOURCES
    comm/Communicator.mpi.cc
    comm/ScopedMpiInit.mpi.cc
  )
  list(APPEND LIBRARY_DEPS MPI::MPI_CXX)
else()
  list(APPEND SOURCES
    comm/Communicator.serial.cc
    comm/ScopedMpiInit.serial.cc
  )
endif()

add_library(celeritas ${SOURCES})
target_link_libraries(celeritas ${LIBRARY_DEPS})

target_include_directories(celeritas
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS celeritas
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

foreach(_SUBDIR base geom)
  file(GLOB _HEADERS
    "${_SUBDIR}/*.h"
  )
  install(FILES ${_HEADERS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${_SUBDIR}/"
  )
endforeach()

#---------------------------------------------------------------------------##