cmake_minimum_required (VERSION 2.8)

project(CmdStan C CXX)

option(BUILD_CMDSTAN_GUIDE
       "Whether the CmdStan Users Guide should be built"
       ON)
option(BUILD_CMDSTAN_TESTS
       "Whether to build the CmdStan specific unit tests"
       ON)

set(CMAKE_MODULE_PATH "${CmdStan_SOURCE_DIR}/stan/cmake")
set(STAN_PATH "${CmdStan_SOURCE_DIR}/stan")

if(BUILD_CMDSTAN_TESTS)
  set(NEED_GTEST ON)
endif()

if(BUILD_CMDSTAN_GUIDE)
  set(NEED_LATEX ON)
endif()

include(StanDependencies)
add_subdirectory(stan)
include(StanBuildFlags)

set(CmdStan_VERSION_MAJOR ${Stan_VERSION_MAJOR})
set(CmdStan_VERSION_MINOR ${Stan_VERSION_MINOR})
set(CmdStan_VERSION_PATCH ${Stan_VERSION_PATCH})
set(CmdStan_VERSION_INTERFACE "1") # Increment this if creating a new release of CmdStan without a new release of Stan
set(CmdStan_VERSION "${CmdStan_VERSION_MAJOR}.${CmdStan_VERSION_MINOR}.${CmdStan_VERSION_PATCH}-${CmdStan_VERSION_INTERFACE}")


include(CTest)

if(BUILD_CMDSTAN_GUIDE)
  add_subdirectory( src/docs/cmdstan-guide )
endif()

# Top level include
include_directories("${Stan_SOURCE_DIR}/src" "${Stan_BINARY_DIR}/src")
include_directories("${CmdStan_SOURCE_DIR}/src" "${CmdStan_BINARY_DIR}/src")

# print 
set(PRINT_SRCS ${PROJECT_SOURCE_DIR}/src/cmdstan/print.cpp)
add_executable(cmdstan-print ${PRINT_SRCS})
set_target_properties(cmdstan-print PROPERTIES OUTPUT_NAME print)
target_link_libraries(cmdstan-print stan)
set(PRINT_BIN ${CMAKE_CURRENT_BINARY_DIR}/print)

# stanc
set(STANC_SRCS ${PROJECT_SOURCE_DIR}/src/cmdstan/stanc.cpp)
add_executable(cmdstan-stanc ${STANC_SRCS})
set_target_properties(cmdstan-stanc PROPERTIES OUTPUT_NAME stanc)
target_link_libraries(cmdstan-stanc stanc)
set(CMDSTANC_BIN ${CMAKE_CURRENT_BINARY_DIR}/stanc)
set(CMDSTAN_MAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/cmdstan/main.cpp")

# unit tests
if(BUILD_CMDSTAN_TESTS)
  add_subdirectory(src/test)
endif()

