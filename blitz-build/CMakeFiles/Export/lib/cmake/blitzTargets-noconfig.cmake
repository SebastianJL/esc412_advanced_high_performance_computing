#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "blitz" for configuration ""
set_property(TARGET blitz APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(blitz PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libblitz.so.1.0.2"
  IMPORTED_SONAME_NOCONFIG "libblitz.so.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS blitz )
list(APPEND _IMPORT_CHECK_FILES_FOR_blitz "${_IMPORT_PREFIX}/lib/libblitz.so.1.0.2" )

# Import target "blitz-static" for configuration ""
set_property(TARGET blitz-static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(blitz-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libblitz.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS blitz-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_blitz-static "${_IMPORT_PREFIX}/lib/libblitz.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
