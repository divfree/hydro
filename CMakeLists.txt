#cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# Default revision tag
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tools/revision.cpp"
     DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/source")

# Static build flag
if (DEFINED STATIC)
  set(STATIC_BUILD ${STATIC})   
else()
  set(STATIC_BUILD OFF)
endif()

# MPI flag
if (DEFINED MPI)
  set(MPI_ENABLE ${MPI})
else()
  set(MPI_ENABLE OFF)
endif()

# Default module selection
if (NOT DEFINED MODULE_HYDRO_2D)
  set(MODULE_HYDRO_2D ON)
endif()
if (NOT DEFINED MODULE_HYDRO_3D)
  set(MODULE_HYDRO_3D OFF)
endif()

# Module selection definitions
if (MODULE_HYDRO_2D)
  add_definitions("-DMODULE_HYDRO_2D")
endif()
if (MODULE_HYDRO_3D)
  add_definitions("-DMODULE_HYDRO_3D")
endif()
if (MODULE_TEST)
  add_definitions("-DMODULE_TEST")
endif()
if (MODULE_HEAT_STORAGE)
  add_definitions("-DMODULE_HEAT_STORAGE")
endif()

# Static build for MSVC
set(CMAKE_USER_MAKE_RULES_OVERRIDE
   ${CMAKE_CURRENT_SOURCE_DIR}/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_CURRENT_SOURCE_DIR}/cxx_flag_overrides.cmake)

project("hydro")

set(BINARY_NAME "hydro")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  

# OpenMP support                              
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              

# GCC warnings and optimizations
if (UNIX)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -std=c++11")
  endif()
  if (DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")  # remove symbolic information
  endif()
endif()

# Source folders
aux_source_directory(source/common SOURCES_COMMON)
aux_source_directory(source/control SOURCES_CONTROL)
aux_source_directory(source/hydro2dmpi SOURCES_HYDRO_MODULE)
aux_source_directory(source/test_module SOURCES_TEST_MODULE)
aux_source_directory(source/heat_storage_module SOURCES_HEAT_STORAGE_MODULE)
aux_source_directory(source SOURCES_ROOT)

set(SOURCE_FILES
  ${SOURCES_COMMON}
  ${SOURCES_CONTROL}
  ${SOURCES_HYDRO_MODULE}
  ${SOURCES_TEST_MODULE}
  ${SOURCES_HEAT_STORAGE_MODULE}
  ${SOURCES_ROOT}
)

# Add target
add_executable(${BINARY_NAME} ${SOURCE_FILES})
     
# Revision tag target
add_custom_target(REVISION_TAG SOURCES "source/revision.cpp")
add_dependencies(${BINARY_NAME} REVISION_TAG)
if (UNIX)
  add_custom_command(
    TARGET REVISION_TAG
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/update_git_rev.sh
    COMMENT " Update Git revision"
  )
else()
  add_custom_command(
    TARGET REVISION_TAG
    COMMAND rem
    COMMENT " Update Git revision (not implemented for Windows)"
  )
endif()

# C++11
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# Static build
if (STATIC_BUILD)
  SET(BUILD_SHARED_LIBRARIES OFF)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_END_STATIC 1)
  if (UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  else()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
endif()

# Tecplot
#if (DEFINED TECPLOT)
#  set(TECPLOT_DIR ${TECPLOT})
#else()
#  set(TECPLOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/tecio")
#endif()
#include_directories(${TECPLOT_DIR}/include)
#find_library(TECIO_LIBRARY 
#             NAMES tecio
#             PATHS ${TECPLOT_DIR}/lib
#             NO_DEFAULT_PATH
#)
#set(LIBRARIES ${LIBRARIES} ${TECIO_LIBRARY})

# Threads
find_package(Threads REQUIRED)
set(LIBRARIES ${LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# MPI
if (MPI_ENABLE) 
  find_package(MPI REQUIRED)
  
  add_definitions("-DMPI_ENABLE")
  
  set(CMAKE_CXX_COMPILER mpicxx)
  set(CMAKE_C_COMPILER mpicc)
  
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(${BINARY_NAME} ${MPI_LIBRARIES})

  if (MPI_COMPILE_FLAGS)
    set_target_properties(${BINARY_NAME} PROPERTIES
                          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if (MPI_LINK_FLAGS)
    set_target_properties(${BINARY_NAME} PROPERTIES
                          LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

target_link_libraries(${BINARY_NAME} ${LIBRARIES})

install(TARGETS ${BINARY_NAME} RUNTIME DESTINATION $ENV{HOME}/bin)
