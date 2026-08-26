cmake_minimum_required(VERSION 3.15)
project(PhaseFieldFerro 
    VERSION 1.0.0 
    LANGUAGES CXX
)

# Supprime le Warning CMP0144 pour MIMALLOC_ROOT
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type de build par defaut : Release" FORCE)
endif()

# ==============================================================================
# GESTION CENTRALISEE DES DEPENDANCES
# ==============================================================================
set(DEPS_ROOT "$ENV{HOME}/dependances_pour_serveur" 
    CACHE PATH "Racine absolue du dossier des dependances")

set(TRACY_DIR "${DEPS_ROOT}/src/tracy-0.13.1-src")
set(GMSH_ROOT "${DEPS_ROOT}/install/gmsh-sdk")
set(SUITESPARSE_ROOT "${DEPS_ROOT}/install/suitesparse-install")
set(MIMALLOC_ROOT "${DEPS_ROOT}/install/mimalloc-install")
set(EIGEN_ROOT "${DEPS_ROOT}/install/eigen/include")
set(PETSC_ROOT "${DEPS_ROOT}/install/petsc")

# ==============================================================================
# OPTIONS DE COMPILATION
# ==============================================================================
option(ENABLE_PROFILING "Activer le profilage avec Tracy" OFF)
option(USE_PARDISO "Utiliser Intel MKL PARDISO pour les solveurs symetriques" OFF)
option(USE_PETSC "Utiliser PETSc/Hypre (BoomerAMG) pour les solveurs symetriques" OFF)

# ==============================================================================
# RECHERCHE DES LIBRAIRIES EXTERNES
# ==============================================================================
find_package(OpenMP REQUIRED)

if(NOT EXISTS "${EIGEN_ROOT}/Eigen/Core")
    message(FATAL_ERROR "Eigen introuvable dans ${EIGEN_ROOT}")
endif()

find_library(GMSH_LIBRARY NAMES gmsh HINTS "${GMSH_ROOT}/lib" REQUIRED)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h HINTS "${GMSH_ROOT}/include" REQUIRED)

# Recherche manuelle et robuste de SuiteSparse (CHOLMOD)
find_path(CHOLMOD_INCLUDE_DIR NAMES cholmod.h 
    HINTS "${SUITESPARSE_ROOT}/include" "${SUITESPARSE_ROOT}/include/suitesparse" 
    REQUIRED)

find_library(CHOLMOD_LIBRARY NAMES cholmod HINTS "${SUITESPARSE_ROOT}/lib" REQUIRED)
find_library(AMD_LIBRARY NAMES amd HINTS "${SUITESPARSE_ROOT}/lib" REQUIRED)
find_library(COLAMD_LIBRARY NAMES colamd HINTS "${SUITESPARSE_ROOT}/lib" REQUIRED)
find_library(SUITESPARSECONFIG_LIBRARY NAMES suitesparseconfig HINTS "${SUITESPARSE_ROOT}/lib" REQUIRED)

# Creation d'une cible virtuelle propre pour CHOLMOD
if(NOT TARGET PkgConfig::CHOLMOD)
    add_library(PkgConfig::CHOLMOD INTERFACE IMPORTED)
    set_target_properties(PkgConfig::CHOLMOD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CHOLMOD_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${CHOLMOD_LIBRARY};${AMD_LIBRARY};${COLAMD_LIBRARY};${SUITESPARSECONFIG_LIBRARY}"
    )
endif()

find_package(mimalloc CONFIG REQUIRED HINTS "${MIMALLOC_ROOT}/lib/cmake/mimalloc-2.4" "${MIMALLOC_ROOT}/lib/cmake/mimalloc-2.1")

# ==============================================================================
# SOURCES DU CORE (Sans main.cpp, qui ira dans l'exécutable)
# ==============================================================================
set(CORE_SOURCES
    core/Core/src/SimulationManager.cpp
    core/Core/src/SimulationApp.cpp
    core/Core/src/TimeManager.cpp

    core/IO/src/IOManager.cpp
    core/IO/src/Datafile.cpp
    core/IO/src/ConfigTypes.cpp
    core/IO/src/ResultsExporter.cpp
    core/IO/src/internal/Diagnostics.cpp

    core/Mesh/src/Mesh.cpp
    core/Mesh/src/MeshGenerators.cpp

    core/Physics/Core/BoundaryManager.cpp
    core/Physics/src/Core/PhysicsManager.cpp
    core/Physics/src/Core/SystemMetrics.cpp
)

set(TEST_SOURCES
    # --- Integration & Safety ---
    core/Tests/src/test_integration_safety.cpp
    core/Tests/src/test_validation_level3_domain_wall.cpp
    
    # --- Core & Utils ---
    core/Core/tests/test_simulation_manager.cpp
    core/Utils/tests/test_rotation2d.cpp
    
    # --- IO ---
    core/IO/tests/test_datafile.cpp
    
    # --- Mesh ---
    core/Mesh/tests/test_mesh.cpp
    core/Mesh/tests/test_gmsh_reader.cpp
    core/Mesh/tests/test_fem.cpp
    core/Mesh/tests/test_bc.cpp
    
    # --- Materials & Microstructure ---
    core/Materials/tests/test_voronoi.cpp
    core/Materials/tests/test_oriented_material.cpp
    core/Materials/tests/test_material_formulas.cpp
    core/Materials/tests/test_material_factory.cpp
    core/Materials/tests/test_materials_manager.cpp
    core/Materials/tests/test_paper_physics.cpp
    core/Materials/tests/test_validation_level0_derivatives.cpp
    core/Materials/tests/test_validation_level1_local_dynamics.cpp
    core/Materials/tests/test_validation_level2_fracture_coupling.cpp
    
    # --- Physics & Equations ---
    core/Physics/tests/test_system_assembler.cpp
    core/Physics/tests/test_polarization_equation.cpp
    core/Physics/tests/test_mechanics_equation.cpp
    core/Physics/tests/test_electrostatics_equation.cpp
    core/Physics/tests/test_fracture_equation.cpp
)

# ==============================================================================
# CONFIGURATION DE LA LIBRAIRIE CORE
# ==============================================================================
if(ENABLE_PROFILING)
    list(APPEND CORE_SOURCES "${TRACY_DIR}/public/TracyClient.cpp")
    set_source_files_properties("${TRACY_DIR}/public/TracyClient.cpp" PROPERTIES COMPILE_OPTIONS "-w")
endif()

if(USE_PETSC)
    list(APPEND CORE_SOURCES core/Utils/src/PetscKspWrapper.cpp)
endif()

add_library(CoreLib STATIC ${CORE_SOURCES})
set_target_properties(CoreLib PROPERTIES LINKER_LANGUAGE CXX)

set(FASTMATH_SENSITIVE_SOURCES
    core/Physics/src/Mechanics/MechanicsAssembler.cpp
    core/Physics/src/Polarization/PolarizationAssembler.cpp
    core/Physics/src/Fracture/FractureAssembler.cpp
    core/Physics/src/Electrostatics/ElectrostaticsAssembler.cpp
)
set_source_files_properties(${FASTMATH_SENSITIVE_SOURCES} PROPERTIES COMPILE_OPTIONS "-fno-fast-math")

target_include_directories(CoreLib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${GMSH_INCLUDE_DIR}
    ${CHOLMOD_INCLUDE_DIR}
    ${EIGEN_ROOT}
)

target_include_directories(CoreLib SYSTEM PUBLIC
    "${TRACY_DIR}/public"
)

target_link_directories(CoreLib PUBLIC "${SUITESPARSE_ROOT}/lib")

target_compile_definitions(CoreLib PUBLIC 
    EIGEN_USE_CHOLMOD 
    EIGEN_NO_DEBUG
)

if(USE_PARDISO)
    if(NOT DEFINED ENV{MKLROOT})
        message(FATAL_ERROR "USE_PARDISO=ON mais MKLROOT n'est pas defini. (Faire 'module load mkl/latest' avant)")
    endif()
    set(MKLROOT "$ENV{MKLROOT}")
    message(STATUS "PARDISO active - MKLROOT = ${MKLROOT}")

    target_compile_definitions(CoreLib PUBLIC USE_PARDISO)
    target_include_directories(CoreLib PUBLIC "${MKLROOT}/include")
    target_link_directories(CoreLib PUBLIC "${MKLROOT}/lib/intel64")

    target_link_libraries(CoreLib PUBLIC
        mkl_intel_lp64
        mkl_gnu_thread
        mkl_core
        gomp
        pthread
        m
        dl
    )
endif()

if(USE_PETSC)
    if(NOT EXISTS "${PETSC_ROOT}/include/petscksp.h")
        message(FATAL_ERROR "USE_PETSC=ON mais PETSc introuvable dans ${PETSC_ROOT}. Lancer scripts/build_petsc.sh d'abord.")
    endif()
    message(STATUS "PETSc/Hypre active - PETSC_ROOT = ${PETSC_ROOT}")

    find_library(PETSC_LIBRARY NAMES petsc HINTS "${PETSC_ROOT}/lib" REQUIRED)

    # ---> La détection officielle de CMake pour MPI <---
    find_package(MPI REQUIRED)

    target_compile_definitions(CoreLib PUBLIC USE_PETSC)
    target_include_directories(CoreLib PUBLIC "${PETSC_ROOT}/include")
    target_link_directories(CoreLib PUBLIC "${PETSC_ROOT}/lib")
    
    # On lie PETSc ET la librairie MPI C++
    target_link_libraries(CoreLib PUBLIC 
        ${PETSC_LIBRARY} 
        MPI::MPI_CXX
    )
endif()

if(ENABLE_PROFILING)
    target_compile_definitions(CoreLib PUBLIC TRACY_ENABLE)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(CoreLib PUBLIC 
    mimalloc
    ${GMSH_LIBRARY}
    PkgConfig::CHOLMOD
) 

target_compile_options(CoreLib PUBLIC -Wall -Wextra -Wpedantic)
target_compile_options(CoreLib PUBLIC
    $<$<CONFIG:Debug>:-O0 -g -fno-omit-frame-pointer>
    $<$<CONFIG:Release>:-O3 -march=native -ffast-math>
    $<$<CONFIG:RelWithDebInfo>:-O3 -march=native -ffast-math -g>
)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
set(ipo_enabled_configs "Release" "RelWithDebInfo")
if(ipo_supported AND CMAKE_BUILD_TYPE IN_LIST ipo_enabled_configs)
    set_target_properties(CoreLib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ==============================================================================
# CREATION DE L'EXECUTABLE (MAIN + TESTS)
# ==============================================================================
# L'exécutable compile main.cpp et la liste complète des tests, puis link la librairie CoreLib
add_executable(PhaseFieldFerro core/Core/src/main.cpp ${TEST_SOURCES})

target_link_libraries(PhaseFieldFerro PRIVATE CoreLib)
set_target_properties(PhaseFieldFerro PROPERTIES
    BUILD_RPATH "${GMSH_ROOT}/lib;${SUITESPARSE_ROOT}/lib;${MIMALLOC_ROOT}/lib;${PETSC_ROOT}/lib"
    INSTALL_RPATH "${GMSH_ROOT}/lib;${SUITESPARSE_ROOT}/lib;${MIMALLOC_ROOT}/lib;${PETSC_ROOT}/lib"
)

if(ipo_supported AND CMAKE_BUILD_TYPE IN_LIST ipo_enabled_configs)
    set_target_properties(PhaseFieldFerro PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_custom_target(clean-output
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMENT "Nettoyage des dossiers de resultats termine !"
)
add_custom_target(clean-export
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean-output
    COMMENT "Alias de clean-output"
)