#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Interpolation::Interpolation" for configuration "Release"
set_property(TARGET Interpolation::Interpolation APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Interpolation::Interpolation PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libInterpolation.so"
  IMPORTED_SONAME_RELEASE "libInterpolation.so"
  )

list(APPEND _cmake_import_check_targets Interpolation::Interpolation )
list(APPEND _cmake_import_check_files_for_Interpolation::Interpolation "${_IMPORT_PREFIX}/lib/libInterpolation.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
