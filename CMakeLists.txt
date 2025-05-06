cmake_minimum_required(VERSION 3.16)
project(SprueCpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find CGAL and its dependencies
# Ensure CGAL is installed and findable by CMake (e.g., via environment variables or CMake prefix path)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Boost REQUIRED COMPONENTS program_options thread system) # For command line and CGAL
find_package(Eigen3 REQUIRED NO_MODULE) # Often needed by CGAL components

# Check if GMP and MPFR are needed (often required for robust kernels)
# find_package(GMP REQUIRED)
# find_package(MPFR REQUIRED)

# Add include directories
include_directories(
    ${CGAL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    # ${GMP_INCLUDE_DIRS}  # Uncomment if GMP found and needed
    # ${MPFR_INCLUDE_DIRS} # Uncomment if MPFR found and needed
    src # Include our source directory
)

# Add executable
add_executable(sprue src/main.cpp)

# Link libraries
target_link_libraries(sprue
    PRIVATE
    CGAL::CGAL # Core CGAL library
    Boost::program_options # For argument parsing
    Boost::thread
    Boost::system
    Eigen3::Eigen
    # ${GMP_LIBRARIES}   # Uncomment if GMP found and needed
    # ${MPFR_LIBRARIES}  # Uncomment if MPFR found and needed
)

# Optional: Enable C++17 features if not default
# target_compile_features(sprue_cpp PRIVATE cxx_std_17)

# Ensure header files from src are accessible if needed later
target_include_directories(sprue PUBLIC src)

# --- Installation (Optional) ---
# install(TARGETS sprue_cpp DESTINATION bin)