#pragma once

#include "Core/Mesh.h"
#include <string>

class MeshGeneratorGmsh {
public:
    static Mesh load_from_msh(const std::string& msh_path, double Lx, double Ly);
};