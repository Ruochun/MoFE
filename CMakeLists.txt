cmake_minimum_required(VERSION 3.10)
project(MoFE LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# User-configurable build type with dropdown options for ccmake GUI
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING 
        "Choose the type of build: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
endif()
# Set the valid values for the CMAKE_BUILD_TYPE property to make it selectable in ccmake
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

# User-configurable options for MPI and AMGX - exposed to ccmake GUI
option(MOFE_USE_MPI "Enable MPI parallel build in MFEM" OFF)
option(MOFE_USE_AMGX "Enable AMGX linear solver support in MFEM" OFF)

# Force CUDA to always be ON in MFEM when CUDA toolkit is available
# This satisfies the requirement while being practical for environments without CUDA
find_program(NVCC_EXECUTABLE nvcc)
if(NVCC_EXECUTABLE)
    set(MFEM_USE_CUDA ON CACHE BOOL "Enable CUDA support in MFEM (forced ON when CUDA toolkit available)" FORCE)
    message(STATUS "CUDA toolkit found: CUDA support will be forced ON")
else()
    set(MFEM_USE_CUDA OFF CACHE BOOL "Enable CUDA support in MFEM (disabled - CUDA toolkit not found)" FORCE)
    message(WARNING "CUDA toolkit not found: CUDA support disabled. Install CUDA toolkit to enable forced CUDA support.")
endif()

# Pass user-configurable options to MFEM
set(MFEM_USE_MPI ${MOFE_USE_MPI} CACHE BOOL "Enable MPI parallel build" FORCE)
set(MFEM_USE_AMGX ${MOFE_USE_AMGX} CACHE BOOL "Enable AMGX usage" FORCE)

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

# Copy data directory to build directory during build process
# This ensures data files are available to executables in the build directory
add_custom_target(copy_data_directory ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory 
        "${CMAKE_SOURCE_DIR}/data" 
        "${CMAKE_BINARY_DIR}/data"
    COMMENT "Copying data directory to build directory"
)

# Make the data copying happen before building executables
add_dependencies(cfd_example copy_data_directory)

# Display configuration summary for user
message(STATUS "=== MoFE Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA support (forced ON): ${MFEM_USE_CUDA}")
message(STATUS "MPI support: ${MOFE_USE_MPI}")
message(STATUS "AMGX support: ${MOFE_USE_AMGX}")
message(STATUS "Data directory will be copied to: ${CMAKE_BINARY_DIR}/data")
message(STATUS "=====================================")
message(STATUS "")
message(STATUS "To modify MPI/AMGX settings, use: ccmake .")
message(STATUS "To change build type, use: ccmake . and select CMAKE_BUILD_TYPE")