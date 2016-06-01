cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Static build flag
if (DEFINED STATIC)
  set(STATIC_BUILD ${STATIC})   
else()
  set(STATIC_BUILD ON)
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

# Boost libraries
set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED 
  COMPONENTS system
             filesystem
             timer
             thread
             chrono
             date_time 
             atomic
)
link_directories(
  ${Boost_LIBRARY_DIRS}
)
include_directories(
  ${Boost_INCLUDE_DIRS}
)
                     
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")
endif()

# Source folders
aux_source_directory(source/common SOURCES_COMMON)
aux_source_directory(source/control SOURCES_CONTROL)
aux_source_directory(source/hydro2dmpi SOURCES_HYDRO_MODULE)
aux_source_directory(source/test_module SOURCES_TEST_MODULE)
aux_source_directory(source SOURCES_ROOT)

set(SOURCE_FILES
  ${SOURCES_COMMON}
  ${SOURCES_CONTROL}
  ${SOURCES_HYDRO_MODULE}
  ${SOURCES_TEST_MODULE}
  ${SOURCES_ROOT}
)

# Add target
add_executable(${BINARY_NAME} ${SOURCE_FILES})

# C++11
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# Static build
if (STATIC_BUILD)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_END_STATIC 1)
endif()

# Tecplot
if (DEFINED TECPLOT)
  set(TECPLOT_DIR ${TECPLOT})
else()
  set(TECPLOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/tecio")
endif()
include_directories(
  ${TECPLOT_DIR}/include
)
find_library(TECIO_LIBRARY 
             NAMES tecio.lib tecio.a
             PATHS ${TECPLOT_DIR}/lib
             NO_DEFAULT_PATH)

find_package(Threads)

# Revision tag
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

set(LIBRARIES
  ${CMAKE_THREAD_LIBS_INIT}
  ${Boost_LIBRARIES}
  ${TECIO_LIBRARY}
)

target_link_libraries(${BINARY_NAME} ${LIBRARIES}) 
