#pragma once
#include "Core/Mesh.h"
#include "IO/Datafile.h" 

class MeshGenerator {
public:
    static Mesh generate_structured_mesh(const MeshConfig& config);
};