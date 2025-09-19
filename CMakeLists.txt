cmake_minimum_required(VERSION 3.10)
project(MoFE LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add debug information
set(CMAKE_BUILD_TYPE Debug)

# Add MFEM subdirectory
add_subdirectory(mfem)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/mfem)
include_directories(${CMAKE_SOURCE_DIR}/mfem)

# Create include and src directories
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/examples)

# Add source files
set(MOFE_SOURCES
    src/cfd_solver.cpp
)

set(MOFE_HEADERS
    include/cfd_solver.hpp
)

# Create MoFE library
add_library(mofe ${MOFE_SOURCES} ${MOFE_HEADERS})
target_link_libraries(mofe mfem)

# Add example executable
add_executable(cfd_example examples/cfd_example.cpp)
target_link_libraries(cfd_example mofe mfem)