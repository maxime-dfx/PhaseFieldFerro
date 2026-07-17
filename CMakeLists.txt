# ==============================================================================
# Configuration Initiale
# ==============================================================================
cmake_minimum_required(VERSION 3.15)
project(PhaseFieldFerro 
    VERSION 1.0.0 
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type de build par defaut : Release" FORCE)
endif()

# ==============================================================================
# Dépendances 
# ==============================================================================
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED) 

# Ajout de Gmsh (Indispensable pour la génération via API)
find_library(GMSH_LIBRARY NAMES gmsh REQUIRED)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h REQUIRED)

# ==============================================================================
# Architecture : Bibliothèque Centrale (CoreLib)
# ==============================================================================
set(CORE_SOURCES
    core/src/Core/BoundaryManager.cpp 

    core/src/Core/Mesh.cpp               
    core/src/Core/MeshGenerator.cpp
    core/src/Core/MeshGeneratorGmsh.cpp

    core/src/IO/Datafile.cpp
    core/src/IO/Diagnostics.cpp
    core/src/IO/MetricsLogger.cpp
    core/src/IO/ResultsExporter.cpp

    core/src/Physics/Fracture.cpp
    core/src/Physics/Polarization.cpp
    core/src/Physics/Mechanics.cpp
    core/src/Physics/Electrostatics.cpp
    core/src/Physics/FractureAssembler.cpp
    core/src/Physics/PolarizationAssembler.cpp
    core/src/Physics/MechanicsAssembler.cpp
    core/src/Physics/ElectrostaticsAssembler.cpp

    core/src/Solvers/LinearSolver.cpp

    core/src/Utils/ProgressBar.cpp

    core/src/Simulation.cpp
    core/src/Validation.cpp
)

add_library(CoreLib STATIC ${CORE_SOURCES})

# ==============================================================================
# Paramétrage de la bibliothèque
# ==============================================================================
set_target_properties(CoreLib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(CoreLib PUBLIC 
    core/include
    ${GMSH_INCLUDE_DIR}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()

# On lie Eigen3 et Gmsh
target_link_libraries(CoreLib PUBLIC 
    Eigen3::Eigen
    ${GMSH_LIBRARY}
) 

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(CoreLib PUBLIC -O3 -march=native -ffast-math -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(CoreLib PUBLIC /O2 /arch:AVX2 /fp:fast /W4)
endif()

# ==============================================================================
# Exécutables
# ==============================================================================

# 1. Application Principale
add_executable(PhaseFieldFerro core/src/main.cpp)
target_link_libraries(PhaseFieldFerro PRIVATE CoreLib)

# 2. Application de Validation
add_executable(Validation core/src/main_validation.cpp)
target_link_libraries(Validation PRIVATE CoreLib)

# ==========================================
# Sécurisation de la mémoire sous Windows (MinGW)
# ==========================================
if(MINGW)
    # Augmente la taille de la stack pour éviter les stack overflows sur les gros maillages
    if(TARGET PhaseFieldFerro)
        target_link_options(PhaseFieldFerro PRIVATE -Wl,--stack,16777216)
    endif()
    if(TARGET Validation)
        target_link_options(Validation PRIVATE -Wl,--stack,16777216)
    endif()
endif()

# ==========================================
# Cibles personnalisées : Nettoyage des datas
# ==========================================
add_custom_target(clean-output
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMENT "Nettoyage des dossiers de resultats termine !"
)

add_custom_target(clean-export
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    # Note : Il n'est pas sécurisé de demander à CMake de supprimer son propre dossier build
    # en cours d'exécution. Si tu veux supprimer build, fais un `rm -rf build` depuis le terminal.
    COMMENT "Nettoyage des outputs termine !"
)