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

# Ajout d'une option pour activer/désactiver le profilage
option(ENABLE_PROFILING "Activer le profilage avec Tracy" OFF)

# ==============================================================================
# Dépendances 
# ==============================================================================
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED) 
find_package(PkgConfig REQUIRED)

# Ajout de Gmsh
find_library(GMSH_LIBRARY NAMES gmsh REQUIRED)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h REQUIRED)

# Ajout de SuiteSparse (CHOLMOD) via pkg-config pour MSYS2
pkg_check_modules(CHOLMOD REQUIRED cholmod)

# ==============================================================================
# Architecture : Bibliothèque Centrale (CoreLib)
# ==============================================================================
set(CORE_SOURCES
    core/src/BC/BoundaryManager.cpp 

    core/src/IO/Datafile.cpp
    core/src/IO/Diagnostics.cpp
    core/src/IO/IOManager.cpp
    core/src/IO/MetricsLogger.cpp
    core/src/IO/ResultsExporter.cpp

    core/src/Mesh/Mesh.cpp               
    core/src/Mesh/MeshGenerator.cpp
    core/src/Mesh/MeshGeneratorGmsh.cpp

    core/src/Physics/Fracture.cpp
    core/src/Physics/Polarization.cpp
    core/src/Physics/Mechanics.cpp
    core/src/Physics/Electrostatics.cpp
    core/src/Physics/FractureAssembler.cpp
    core/src/Physics/PolarizationAssembler.cpp
    core/src/Physics/MechanicsAssembler.cpp
    core/src/Physics/ElectrostaticsAssembler.cpp

    core/src/Simulation/TimeManager.cpp
    core/src/Simulation/Simulation.cpp

    core/src/Solvers/LinearSolver.cpp

    core/src/Utils/ProgressBar.cpp

    core/src/Validation.cpp
)

# Ajout du client Tracy aux sources si le profilage est activé
if(ENABLE_PROFILING)
    list(APPEND CORE_SOURCES "C:/Users/amorvan/Documents/Librairies/tracy/public/TracyClient.cpp")
endif()

add_library(CoreLib STATIC ${CORE_SOURCES})

# ==============================================================================
# Paramétrage et Optimisations Extrêmes de la bibliothèque
# ==============================================================================
set_target_properties(CoreLib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(CoreLib PUBLIC 
    core/include
    ${GMSH_INCLUDE_DIR}
    ${CHOLMOD_INCLUDE_DIRS}
)

if(ENABLE_PROFILING)
    target_include_directories(CoreLib PUBLIC "C:/Users/amorvan/Documents/Librairies/tracy/public")
endif()

# Activation d'Eigen CHOLMOD et désactivation des assertions de sécurité (Vitesse maximale)
target_compile_definitions(CoreLib PUBLIC 
    EIGEN_USE_CHOLMOD 
    EIGEN_NO_DEBUG
)

# Activation du flag de préprocesseur pour Tracy
if(ENABLE_PROFILING)
    target_compile_definitions(CoreLib PUBLIC TRACY_ENABLE)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Liaison Eigen3, Gmsh et CHOLMOD
target_link_libraries(CoreLib PUBLIC 
    Eigen3::Eigen
    ${GMSH_LIBRARY}
    ${CHOLMOD_LIBRARIES}
) 

# Liaison des bibliothèques système Windows requises par Tracy sous MSYS2/MinGW
if(ENABLE_PROFILING AND MINGW)
    target_link_libraries(CoreLib PUBLIC ws2_32 dbghelp psapi)
endif()

# Flags de compilation agressifs
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(CoreLib PUBLIC -O3 -march=native -ffast-math -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(CoreLib PUBLIC /O2 /arch:AVX2 /fp:fast /W4)
endif()

# ==============================================================================
# Optimisation globale LTO/IPO (Link Time Optimization)
# ==============================================================================
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_target_properties(CoreLib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ==============================================================================
# Exécutables
# ==============================================================================

# 1. Application Principale
add_executable(PhaseFieldFerro core/src/main.cpp)
target_link_libraries(PhaseFieldFerro PRIVATE CoreLib)
if(ipo_supported)
    set_target_properties(PhaseFieldFerro PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# 2. Application de Validation
add_executable(Validation core/src/main_validation.cpp)
target_link_libraries(Validation PRIVATE CoreLib)
if(ipo_supported)
    set_target_properties(Validation PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ==========================================
# Sécurisation de la mémoire sous Windows (MSYS2/MinGW)
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
    COMMENT "Nettoyage des outputs termine !"
)