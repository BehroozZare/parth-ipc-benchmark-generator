# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if (HAS_PARENT)
    set(IPC_TOPLEVEL_PROJECT OFF)
else ()
    set(IPC_TOPLEVEL_PROJECT ON)
endif ()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if (IPC_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else ()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if (${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build IPC is ${REQUIRED_CMAKE_VERSION}")
    endif ()
endif ()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/IPCOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/IPCOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/IPCOptions.cmake)
endif ()

################################################################################

project(IPC DESCRIPTION "Incremental Potential Contact (IPC)" LANGUAGES CXX)

################################################################################

# Project options
option(IPC_WITH_OPENGL "Use libgil viewer for preview and output" ON)
option(IPC_WITH_AMGCL "Enable AMGCL linear system solver" OFF)
option(IPC_WITH_CHOLMOD "Enable CHOLMOD linear system solver" ON)

option(IPC_WITH_STRUMPACK "Use STRUMPACK as the solver" OFF)
option(IPC_WITH_PARTH "Use PARTH as the solver" OFF)
option(IPC_WITH_MKL "Enable MKL linear system solver" OFF)
option(IPC_WITH_METIS "Enable METIS for PARTH integration" ON)
option(IPC_WITH_SYMPILER "Enable SYMPILER linear system solver" OFF)

set(IPC_DEFAULT_LINSYSSOLVER "CHOLMOD" CACHE STRING "Linear system solver to use (options: CHOLMOD, AMGCL, EIGEN)")
set_property(CACHE IPC_DEFAULT_LINSYSSOLVER PROPERTY STRINGS "CHOLMOD" "AMGCL" "EIGEN")
## SQP and QP options
option(IPC_WITH_GUROBI "Enable Gurobi QP solver for SQP method" OFF)
option(IPC_WITH_OSQP_MKL "Use MKL Pardiso linear solver in OSQP if available" ON)
option(IPC_EXIT_UPON_QP_FAIL "Exit if the QP/SQP intersects or blows-up" OFF)
## CCD options
option(IPC_WITH_EXACT_CCD "Use exact CCD as a verification" OFF)
## Top level options
option(IPC_WITH_TESTS "Build unit tests using Catch2" ${IPC_TOPLEVEL_PROJECT})


################################################################################

### Configuration
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ipc)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/find)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes)

################################################################################
# IPC Library
################################################################################

# project source files
file(GLOB SRCFILES
        "src/*.cpp"
        "src/Energy/*.cpp"
        "src/Energy/Physics_Elasticity/*.cpp"
        "src/Energy/Geometry_Isometry/*.cpp"
        "src/Energy/Connectivity/*.cpp"
        "src/Energy/Collision/*.cpp"
        "src/CollisionObject/*.cpp"
        "src/CCD/*.cpp"
        "src/Utils/*.cpp"
        "src/LinSysSolver/*.cpp"
        "src/TimeStepper/*.cpp"
        "PARTH/*.cpp"
)

FOREACH (item ${SRCFILES})
    IF (${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF (${item} MATCHES "main.cpp")
ENDFOREACH (item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC
        "src"
        "src/Energy"
        "src/Energy/Physics_Elasticity"
        "src/Energy/Geometry_Isometry"
        "src/Energy/Connectivity"
        "src/Energy/Collision"
        "src/CollisionObject/"
        "src/CCD/"
        "src/Utils"
        "src/Utils/SVD"
        "src/LinSysSolver"
        "src/TimeStepper"
        "PARTH/include"
)

################################################################################
# Required Libraries
################################################################################
if (IPC_WITH_METIS)
    find_package(METIS REQUIRED)
    if (EXISTS "${METIS_INCLUDES}")
        message(STATUS "The METIS has been found")
    else ()
        message(STATUS "METIS not found.")
        set(METIS_LIBRARIES "")
        set(MKL 0)
    endif ()

    target_include_directories(${PROJECT_NAME}_dev PUBLIC "${METIS_INCLUDES}")
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC
            ${METIS_LIBRARIES})
endif ()


#MKL
if (IPC_WITH_MKL)
    ### Find Intel MKL if available
    find_package(MKL OPTIONAL_COMPONENTS)
    if (EXISTS "${MKL_INCLUDE_DIRS}")
        message(STATUS "The MKL has been found")
    else ()
        message(STATUS "Intel MKL not found.")
        set(MKL_LIBRARIES "")
        set(MKL 0)
    endif ()

    target_include_directories(${PROJECT_NAME}_dev PUBLIC "${MKL_INCLUDE_DIRS}")
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC
            ${MKL_LIBRARIES})
    add_compile_options(-m64 -W)
endif ()

#Strumpack
if (IPC_WITH_STRUMPACK)
    message(STATUS "FINDING STRUMPACK")
    find_package(MPI REQUIRED)
    target_link_libraries(${PROJECT_NAME}_dev
            PUBLIC
            MPI::MPI_CXX)
    set(STRUMPACK_DIR "$ENV{STRUMPACKROOT}/lib/cmake/STRUMPACK/")
    find_package(STRUMPACK REQUIRED)

    target_include_directories(${PROJECT_NAME}_dev PUBLIC "$ENV{STRUMPACKROOT}/include")
    target_link_libraries(${PROJECT_NAME}_dev PRIVATE STRUMPACK::strumpack)
endif ()

# SuiteSparse
if (IPC_WITH_CHOLMOD)
    #    if (DEFIf
    set(SUITESPARSE_INCLUDE_DIR_HINTS $ENV{SUITESPARSE_INC})
    set(SUITESPARSE_LIBRARY_DIR_HINTS $ENV{SUITESPARSE_LIB})
    find_package(SuiteSparse REQUIRED)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_LIBRARIES})
    target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_INCLUDE_DIRS})
endif ()



# Extra warnings
include(ipc_warnings)
target_link_libraries(${PROJECT_NAME}_dev PRIVATE IPC::warnings)

# tbb
include(onetbb)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC TBB::tbb)


#filesystem
if (APPLE)
    # On Apple, we don't need to link stdc++fs as filesystem is part of the standard library
    # The GHC filesystem library will handle the implementation
else()
    # On other platforms, link stdc++fs for filesystem support
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC stdc++fs)
endif()

# OSQP
include(osqp)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC osqp::osqp)

## libigl
#include(libigl)
## Other modules you could enable
##igl_include(embree)
#igl_include(imgui)
#igl_include(stb)
##igl_include(opengl)
##igl_include(stb)
##igl_include(predicates)
##igl_include(xml)
##igl_include(copyleft cgal)
##igl_include(copyleft comiso)
#igl_include(copyleft core)
##igl_include(copyleft cork)
#igl_include(copyleft tetgen)
##igl_include(restricted matlab)
##igl_include(restricted mosek)
#igl_include(restricted triangle)
#
#target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl_copyleft::core igl_restricted::triangle igl_copyleft::tetgen)
#if (IPC_WITH_OPENGL)
#    igl_include(glfw)
#    target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::glfw igl::imgui igl::stb)
#    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_OPENGL)
#endif ()

# libigl
include(eigen)
include(libigl)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::core)
if (IPC_WITH_OPENGL)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::glfw igl::imgui igl::stb)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_OPENGL)
endif ()


if (LIBIGL_WITH_PREDICATES)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC igl::predicates)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_PREDICATES)
endif ()

# CCD
include(ccd_wrapper)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ccd_wrapper::ccd_wrapper)
if (IPC_WITH_EXACT_CCD)
    message(WARNING "The exact CCD methods RootParity and BSC miss some collisions! Please use TightInclusion instead.")
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_EXACT_CCD)
endif ()
if (CCD_WRAPPER_WITH_FPRP)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_FPRP_CCD)
endif ()
if (CCD_WRAPPER_WITH_TIGHT_INCLUSION)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_TIGHT_INCLUSION)
endif ()

# Logger
include(spdlog)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC spdlog::spdlog)

# AMGCL
if (IPC_WITH_AMGCL)
    include(amgcl)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC amgcl::amgcl)
endif ()

if (IPC_WITH_GUROBI)
    include(eigen_gurobi)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC EigenGurobi::EigenGurobi)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_GUROBI)
endif ()

# MshIO
include(mshio)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC mshio::mshio)

# Filesystem
include(filesystem)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ghc::filesystem)

# target_link_libraries(${PROJECT_NAME}_dev PUBLIC pardiso600-GNU720-X86-64 gfortran)
# find_package(Pardiso)
# if(PARDISO_FOUND)
#     find_package(LAPACK)
#     if(LAPACK_FOUND)
#         target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${LAPACK_LIBRARIES})
#     else()
#         message(FATAL_ERROR “unable to find lapack”)
#     endif()


find_package(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DOPENMP)
endif ()

if (APPLE)
    set(OMP_BREW /usr/local/Cellar/libomp/18.1.4/)
    target_include_directories(${PROJECT_NAME}_dev PUBLIC ${OMP_BREW}/include)
    find_library(OMP NAMES omp HINTS ${OMP_BREW}/lib)
    if (OMP)
        message(STATUS "Found omp: ${OMP}")
    else ()
        message(FATAL_ERROR "omp library not found")
    endif ()
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${OMP})
    add_definitions(-DOPENMP)
endif ()


#     target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${PARDISO_LIBRARIES})
# else()
#     message(WARNING “Pardiso not found, solver will not be available.“)
# endif()
# target_link_libraries(${PROJECT_NAME}_dev PUBLIC gfortran)

################################################################################
# Simulations Executable
################################################################################

add_executable(${PROJECT_NAME}_bin "src/main.cpp")
# Link against the library
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)
# CLI11
include(cli11)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC CLI11::CLI11)

################################################################################
# Sub-Projects
################################################################################

add_subdirectory(src/Projects)

################################################################################
# Code Format Tool
################################################################################

find_program(CLANG_FORMAT_INSTALLED "clang-format")
if (CLANG_FORMAT_INSTALLED)
    message(STATUS "Building code format tool")
    add_subdirectory(Format)
else ()
    message(STATUS "clang-format not found: skipping Format subdirectory")
endif ()


################################################################################
# Compiler options
################################################################################

# Set version of SVD to use
message(STATUS "Using scalar implmentation of SVD")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC -DUSE_SCALAR_IMPLEMENTATION)

# Figure out SSE level support
message(STATUS "Seaching for SSE...")
find_package(SSE)

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)

# Figure out FMA level support
message(STATUS "Searching for FMA...")
find_package(FMA)

# Add SSE, AVX, and FMA flags to compiler flags
string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
target_compile_options(${PROJECT_NAME}_dev PUBLIC ${SIMD_FLAGS})

# Add -pthread to compilation and linking
find_package(Threads)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${CMAKE_THREAD_LIBS_INIT})

# Use C++17
target_compile_features(${PROJECT_NAME}_dev PUBLIC cxx_std_17)

# code only support 3D
message(STATUS "Building for three dimensions")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC DIM=3)

# Ensure that PRId64 is defined
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC -D__STDC_FORMAT_MACROS)

# Unit tests
if (IPC_WITH_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif ()

# Developing
if (IPC_WITH_PARTH)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_PARTH)
endif ()


# MKL Pardiso OSQP linear system solver
if (IPC_WITH_OSQP_MKL)
    find_library(MKL_RT_LIBRARY mkl_rt
            PATHS $ENV{MKLROOT}/lib
            $ENV{MKLROOT}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            /opt/intel/mkl/lib)
    if (MKL_RT_LIBRARY)
        message(STATUS "MKL found: enabling OSQP MKL Pardiso linear system solver")
        target_compile_definitions(${PROJECT_NAME}_dev PUBLIC OSQP_USE_MKL_PARDISO)
    else ()
        message(WARNING "MKL not found: disabling OSQP MKL Pardiso linear system solver")
    endif ()
endif ()

# Set the selected linear system solver
if (IPC_DEFAULT_LINSYSSOLVER STREQUAL "CHOLMOD" AND NOT IPC_WITH_CHOLMOD)
    message(FATAL_ERROR "CHOLMOD selected as IPC_DEFAULT_LINSYSSOLVER, but CHOLMOD is disabled (IPC_WITH_CHOLMOD=OFF)!")
endif ()
if (IPC_DEFAULT_LINSYSSOLVER STREQUAL "AMGCL" AND NOT IPC_WITH_AMGCL)
    message(FATAL_ERROR "AMGCL selected as IPC_DEFAULT_LINSYSSOLVER, but AMGCL is disabled (IPC_WITH_AMGCL=OFF)!")
endif ()
string(TOUPPER "${IPC_DEFAULT_LINSYSSOLVER}" IPC_DEFAULT_LINSYSSOLVER)
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC "IPC_DEFAULT_LINSYSSOLVER=LinSysSolverType::${IPC_DEFAULT_LINSYSSOLVER}")

if (IPC_WITH_AMGCL)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_AMGCL)
endif ()
if (IPC_WITH_CHOLMOD)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_CHOLMOD)
endif ()
if (IPC_WITH_MKL)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_MKL)
    #    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EIGEN_USE_MKL_ALL)
endif ()
if (IPC_WITH_STRUMPACK)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IPC_WITH_STRUMPACK)
endif ()

# Quiet the viewer
if (LIBIGL_QUIET_VIEWER)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC IGL_VIEWER_VIEWER_QUIET)
endif ()

if (IPC_EXIT_UPON_QP_FAIL)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EXIT_UPON_QP_FAIL)
endif ()
