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
    core/src/Physics/Math.cpp
    core/src/Physics/Material.cpp

    core/src/Solvers/LinearSolver.cpp

    core/src/Utils/ProgressBar.cpp

    core/src/Simulation.cpp
    

    # --- Vos nouveaux fichiers Header-Only (Optionnels pour la compilation, utiles pour l'IDE) ---
    # core/tests/FerroelectricMaterial.hpp
    # core/tests/StaggeredSolver.hpp
    # core/tests/BoundaryManager.hpp
    # core/tests/MultiphysicsAssembler.hpp
    # core/tests/Datafile.hpp
    # core/tests/GmshReader.hpp
    # core/tests/Element.hpp
)

add_library(CoreLib STATIC ${CORE_SOURCES})

# ==============================================================================
# FORCER LE LANGAGE C++ (INDISPENSABLE CAR IL N'Y A QUE DES .hpp AU-DESSUS)
# ==============================================================================
set_target_properties(CoreLib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(CoreLib PUBLIC core/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(CoreLib PUBLIC Eigen3::Eigen) 

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

# 2. Unités de Tests
# add_executable(UnitTests core/tests/test_physics.cpp)
# target_link_libraries(UnitTests PRIVATE CoreLib)

# add_executable(MeshTests core/tests/test_mesh.cpp)
# target_link_libraries(MeshTests PRIVATE CoreLib)

# enable_testing()
# add_test(NAME PhysicsTest COMMAND UnitTests)
# add_test(NAME MeshTest COMMAND MeshTests)

# ==========================================
# Sécurisation de la mémoire sous Windows (MinGW)
# ==========================================
if(MINGW)
    # Vérifie si la cible existe avant d'appliquer l'option (évite les crashs CMake)
    if(TARGET PhaseFieldFerro)
        target_link_options(PhaseFieldFerro PRIVATE -Wl,--stack,16777216)
    endif()
    if(TARGET UnitTests)
        target_link_options(UnitTests PRIVATE -Wl,--stack,16777216)
    endif()
    if(TARGET MeshTests)
        target_link_options(MeshTests PRIVATE -Wl,--stack,16777216)
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
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/build"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/build"
    COMMENT "Nettoyage complet (output et build) termine !"
)