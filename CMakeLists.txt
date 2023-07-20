cmake_minimum_required(VERSION 3.14)
project(my_project)

set(CMAKE_CXX_STANDARD 17)
cmake_policy(SET CMP0135 NEW)

# Add the following lines to include Google Test via FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  my_project_test
  test.cpp
  s21_Operations.cpp
  matrix.cc
)

target_link_libraries(
  my_project_test
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(my_project_test)
