cmake_minimum_required (VERSION 2.8)

project(CmdStan C CXX)

add_subdirectory(stan)

set(CmdStan_VERSION_MAJOR ${Stan_VERSION_MAJOR})
set(CmdStan_VERSION_MINOR ${Stan_VERSION_MINOR})
set(CmdStan_VERSION_PATCH ${Stan_VERSION_PATCH})
set(CmdStan_VERSION_INTERFACE "1") # Increment this if creating a new release of CmdStan without a new release of Stan
set(CmdStan_VERSION "${CmdStan_VERSION_MAJOR}.${CmdStan_VERSION_MINOR}.${CmdStan_VERSION_PATCH}-${CmdStan_VERSION_INTERFACE}")

set(CMAKE_MODULE_PATH ${Stan_SOURCE_DIR}/cmake)

include(CTest)

option(BUILD_CMDSTAN_GUIDE
       "Whether the CmdStan Users Guide should be built"
       ON)

if(BUILD_CMDSTAN_GUIDE)
  include(UseLATEX)

  add_subdirectory( src/docs/cmdstan-guide )
endif()

# Find Eigen
if(USE_SYSTEM_EIGEN)
  find_package( Eigen3 REQUIRED )
  include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}")
else()
  include_directories(SYSTEM "${Stan_SOURCE_DIR}/lib/eigen_3.2.0")
endif()

# Find Boost
if(USE_SYSTEM_BOOST)
  find_package( Boost 1.54 REQUIRED )
  include_directories(SYSTEM "${Boost_INCLUDE_DIR}")
else()
  include_directories(SYSTEM "${Stan_SOURCE_DIR}/lib/boost_1.54.0")
endif()

# Setup GTest library
if(USE_SYSTEM_GTEST)
  find_package( GTest REQUIRED )
  set(GTEST_LIBRARIES ${GTEST_BOTH_LIBRARIES})
endif()

# Setup compiler flags
check_cxx_compiler_flag("-pipe" HAVE_PIPE_FLAG)
if(HAVE_PIPE_FLAG)
  add_compile_options( -pipe )
endif(HAVE_PIPE_FLAG)

# FIXME: If any of these warning flags will break compilers, we can check for
# them like is being done with -pipe
add_compile_options( -Wall -Wno-unused-function -Wno-tautological-compare
                     -Wno-c++11-long-long )

# This is used for the syntax-only test compiles
check_cxx_compiler_flag("-fsyntax-only" HAVE_SYNTAXONLY_FLAG)

# Setup boost defines
add_definitions(-DBOOST_RESULT_OF_USE_TR1
                -DBOOST_NO_DECLTYPE
                -DBOOST_DISABLE_ASSERTS)
# Top level include
include_directories("${Stan_SOURCE_DIR}/src" "${Stan_BINARY_DIR}/src")

# print 
set(PRINT_SRCS ${PROJECT_SOURCE_DIR}/src/cmdstan/print.cpp)
add_executable(print ${PRINT_SRCS})
target_link_libraries(print stan)
#set(PRINT_BIN ${CMAKE_CURRENT_BINARY_DIR}/print)

# stanc
set(STANC_SRCS ${PROJECT_SOURCE_DIR}/src/cmdstan/stanc.cpp)
add_executable(cmdstan-stanc ${STANC_SRCS})
set_target_properties(cmdstan-stanc PROPERTIES OUTPUT_NAME stanc)
target_link_libraries(cmdstan-stanc stanc)
#set(STANC_BIN ${CMAKE_CURRENT_BINARY_DIR}/stanc)

# unit tests
#include(CTest)
#add_subdirectory(src/test)

