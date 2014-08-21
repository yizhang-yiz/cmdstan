include(StanTesting)

include_directories(SYSTEM "${GTEST_INCLUDE_DIRS}")

# Disable optimization for building tests
foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
  string(REGEX REPLACE "-O[^0]" "-O0" ${flag_var} "${${flag_var}}")
  string(REGEX REPLACE "/O[^0]" "/O0" ${flag_var} "${${flag_var}}")
endforeach(flag_var)

add_subdirectory(test-models)
add_subdirectory(models)
add_subdirectory(interface)

