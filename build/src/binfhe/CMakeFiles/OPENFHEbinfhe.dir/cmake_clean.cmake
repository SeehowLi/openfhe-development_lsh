file(REMOVE_RECURSE
  "../../lib/libOPENFHEbinfhe.dll"
  "../../lib/libOPENFHEbinfhe.dll.a"
  "../../lib/libOPENFHEbinfhe.dll.manifest"
  "../../lib/libOPENFHEbinfhe.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/OPENFHEbinfhe.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
