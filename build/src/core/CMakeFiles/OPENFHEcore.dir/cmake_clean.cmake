file(REMOVE_RECURSE
  "../../lib/libOPENFHEcore.dll"
  "../../lib/libOPENFHEcore.dll.a"
  "../../lib/libOPENFHEcore.dll.manifest"
  "../../lib/libOPENFHEcore.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/OPENFHEcore.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
