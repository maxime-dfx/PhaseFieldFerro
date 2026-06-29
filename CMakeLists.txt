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

# ==============================================================================
# Architecture : Bibliothèque Centrale (CoreLib)
# ==============================================================================
set(CORE_SOURCES
    core/src/Core/BoundaryManager.cpp
    core/src/Core/Mesh.cpp

    core/src/IO/Datafile.cpp
    core/src/IO/MetricsLogger.cpp
    core/src/IO/ResultsExporter.cpp

    core/src/Utils/ProgressBar.cpp
    core/src/Utils/Fields.cpp

    core/src/Physics/Material.cpp
    core/src/Physics/Math.cpp
    core/src/Physics/Polarization.cpp
    core/src/Physics/Mechanics.cpp
    core/src/Physics/Electrostatics.cpp
    core/src/Physics/Fracture.cpp

    core/src/Simulation.cpp

    core/src/main.cpp
)

add_library(CoreLib STATIC ${CORE_SOURCES})

# PUBLIC permet aux exécutables qui utilisent CoreLib de trouver les headers
target_include_directories(CoreLib PUBLIC core/include)

# Lier OpenMP et Eigen à la bibliothèque
if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(CoreLib PUBLIC Eigen3::Eigen) # <-- NOUVEAU : Lien avec Eigen

# Flags HPC appliqués à la bibliothèque
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

# 2. Unité de Tests
add_executable(UnitTests core/tests/test_physics.cpp)
target_link_libraries(UnitTests PRIVATE CoreLib)
add_executable(MeshTests core/tests/test_mesh.cpp)
target_link_libraries(MeshTests PRIVATE CoreLib)

enable_testing()
add_test(NAME PhysicsTest COMMAND UnitTests)
add_test(NAME MeshTest COMMAND MeshTests)

# ==========================================
# Cible personnalisée : Nettoyage des datas
# ==========================================
add_custom_target(clean-output
    # Note : remove_directory supprime le dossier entier, on ne met pas de *.vtk à la fin.
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMENT "Nettoyage des dossiers de resultats termine !"
)
add_custom_target(clean-export
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/output"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/build"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/build"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/build"
    COMMENT "Nettoyage des dossiers de resultats termine !"
)