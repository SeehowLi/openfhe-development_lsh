file(REMOVE_RECURSE
  "../../lib/libOPENFHEpke.dll"
  "../../lib/libOPENFHEpke.dll.a"
  "../../lib/libOPENFHEpke.dll.manifest"
  "../../lib/libOPENFHEpke.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/OPENFHEpke.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
