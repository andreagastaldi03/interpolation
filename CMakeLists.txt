cmake_minimum_required(VERSION 3.15)

# ============================
# SETUP
# ============================
project(Interpolation CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "No build type specified. Defaulting to Release.")
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

set(Interpolation_Version 0.0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================
# MAIN
# ============================


file(GLOB source_files
   src/*.cc
)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(Interpolation SHARED ${source_files})
add_library(Interpolation::Interpolation ALIAS Interpolation)

target_compile_features(Interpolation PUBLIC cxx_std_20)

target_include_directories(Interpolation
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/src
)

# ============================
# COMPILE FLAGS
# ============================
target_compile_options(Interpolation
  PUBLIC
    -Wall -Wextra          
  PRIVATE
    -fPIC                  
)

target_compile_options(Interpolation
  PUBLIC
    $<$<CONFIG:Release>:-O3 -funroll-loops>
    $<$<CONFIG:RelWithDebInfo>:-O2 -g -funroll-loops -fsanitize=address,undefined,alignment,bounds,null,vptr>
    $<$<CONFIG:Debug>:-O0 -g -funroll-loops>
)

enable_testing()
add_subdirectory(tests)
add_subdirectory(run)

# ============================
# INSTALLATION
# ============================
include(CMakePackageConfigHelpers)

install(TARGETS Interpolation
    EXPORT InterpolationTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hh"
)

install(
    EXPORT InterpolationTargets
    FILE InterpolationTargets.cmake
    NAMESPACE Interpolation::
    DESTINATION lib/cmake/Interpolation
)

# Install config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/InterpolationConfigVersion.cmake"
    DESTINATION lib/cmake/Interpolation
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/InterpolationConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/InterpolationConfig.cmake
    INSTALL_DESTINATION lib/cmake/Oros
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/InterpolationConfigVersion.cmake
    VERSION 1.0.0
    COMPATIBILITY AnyNewerVersion
)

# Install Config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/InterpolationConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/InterpolationConfigVersion.cmake
    DESTINATION lib/cmake/Oros
)



find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
   message(STATUS "Found Doxygen: ${DOXYGEN_EXECUTABLE}")

   set(DOXYGEN_DIR ${CMAKE_SOURCE_DIR}/docs/doxygen)
   set(DOXYFILE ${DOXYGEN_DIR}/Doxyfile)
   set(DOC_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/docs)

   add_custom_target(Interpolation-doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
      WORKING_DIRECTORY ${DOXYGEN_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM
   )
else()
   message(WARNING "Doxygen not found in the system! Cannot generate documentation.")
endif()

