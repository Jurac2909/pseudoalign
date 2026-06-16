cmake_minimum_required(VERSION 3.16)
project(pseudokall CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O2 -Wall -Wextra)

# Main library (shared between tool and tests)
add_library(pklib
    src/fasta.cpp
    src/dbg.cpp
)
target_include_directories(pklib PUBLIC src)

# Test driver
add_executable(test_dbg tests/test_dbg.cpp)
target_link_libraries(test_dbg pklib)
