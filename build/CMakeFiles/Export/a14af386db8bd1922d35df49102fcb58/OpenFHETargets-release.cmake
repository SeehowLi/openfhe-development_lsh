#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "OPENFHEcore" for configuration "Release"
set_property(TARGET OPENFHEcore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(OPENFHEcore PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEcore.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEcore.dll"
  )

list(APPEND _cmake_import_check_targets OPENFHEcore )
list(APPEND _cmake_import_check_files_for_OPENFHEcore "${_IMPORT_PREFIX}/lib/libOPENFHEcore.dll.a" "${_IMPORT_PREFIX}/lib/libOPENFHEcore.dll" )

# Import target "OPENFHEpke" for configuration "Release"
set_property(TARGET OPENFHEpke APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(OPENFHEpke PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEpke.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEpke.dll"
  )

list(APPEND _cmake_import_check_targets OPENFHEpke )
list(APPEND _cmake_import_check_files_for_OPENFHEpke "${_IMPORT_PREFIX}/lib/libOPENFHEpke.dll.a" "${_IMPORT_PREFIX}/lib/libOPENFHEpke.dll" )

# Import target "OPENFHEbinfhe" for configuration "Release"
set_property(TARGET OPENFHEbinfhe APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(OPENFHEbinfhe PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEbinfhe.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libOPENFHEbinfhe.dll"
  )

list(APPEND _cmake_import_check_targets OPENFHEbinfhe )
list(APPEND _cmake_import_check_files_for_OPENFHEbinfhe "${_IMPORT_PREFIX}/lib/libOPENFHEbinfhe.dll.a" "${_IMPORT_PREFIX}/lib/libOPENFHEbinfhe.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
