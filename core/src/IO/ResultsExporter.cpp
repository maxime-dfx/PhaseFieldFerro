#include <fstream>
#include <iostream>
#include <cstdint>
#include <filesystem>
#include "IO/ResultsExporter.h"
#include "Utils/Logger.h"
#include "Utils/EnsureDir.h" 
#include <Eigen/Dense>


ResultsExporter::ResultsExporter(const Mesh& msh) : m_mesh(msh) {}

namespace fs = std::filesystem;

void ResultsExporter::exportToMesh(const std::string& filename) {
    Logger::info("Demarrage de l'export du maillage vers : " + filename);

    namespace fs = std::filesystem;
    fs::path p(filename);

    EnsureDir ensureDir;
    if (p.has_parent_path()) {
        ensureDir.prepareOutputDirectory(p.parent_path().string());
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // 3. Écriture des données
    out << "MeshVersionFormatted 1\nDimension 3\n\n";
    
    auto nodes = m_mesh.get_nodes();
    out << "Vertices\n" << nodes.size() << "\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << " "  << node.ref_tag << "\n";
    }
    out << "\n";

    if (m_mesh.get_element_type() == ElementType::QUAD4) {
        auto quads = m_mesh.get_elements();
        out << "Quadrilaterals\n" << quads.size() << "\n";
        for (const auto& q : quads) {
            out << q.node_indices[0] << " " << q.node_indices[1] << " " 
                << q.node_indices[2] << " " << q.node_indices[3] << " " << q.ref_tag << "\n";
        }
    } else if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        auto tris = m_mesh.get_elements();
        out << "Triangles\n" << tris.size() << "\n";
        for (const auto& t : tris) {
            out << (t.node_indices[0] + 1) << " " << (t.node_indices[1] + 1) << " " 
                << (t.node_indices[2] + 1) << " " << t.ref_tag << "\n";        
        }
    } else {
        Logger::error("Type d'élément non pris en charge.");
    }
    
    out << "\nEnd\n";
    out.close();
}

void ResultsExporter::exportScalarVTK(const std::string& filename, 
                                      const std::string& field_name,   
                                      const Eigen::VectorXd& data_field) { 

    namespace fs = std::filesystem;
    fs::path p(filename);

    EnsureDir ensureDir;
    if (p.has_parent_path()) {
        ensureDir.prepareOutputDirectory(p.parent_path().string());
    }

    std::ofstream out(filename);
    out.imbue(std::locale("C")); 

    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Simulation\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    auto nodes = m_mesh.get_nodes();
    out << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << "\n";
    }

    if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        auto elements = m_mesh.get_elements();
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 4) << "\n";
        for (const auto& tri : elements) {
            out << "3 " << tri.node_indices[0] << " " << tri.node_indices[1] << " " << tri.node_indices[2] << "\n";
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "5\n"; 
        }
    } 
    else if (m_mesh.get_element_type() == ElementType::QUAD4) {
        auto elements = m_mesh.get_elements();
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 5) << "\n";
        for (const auto& quad : elements) {
            out << "4 " << quad.node_indices[0] << " " << quad.node_indices[1] << " " << quad.node_indices[2] << " " << quad.node_indices[3] << "\n";
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "9\n"; 
        }
    } 
    else {
        Logger::error("Type d'élément non pris en charge pour l'export VTK.");
        return;
    }

    out << "\nPOINT_DATA " << nodes.size() << "\n";
    
    if (data_field.size() != nodes.size()) {
        Logger::error("La taille du champ scalaire '" + field_name + "' ne correspond pas au nombre de noeuds.");
        return; 
    }

    out << "SCALARS " << field_name << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& val : data_field) {
        out << val << "\n";
    }
    
    out.close();
}

void ResultsExporter::exportVectorVTK(const std::string& filename, 
                                                  const std::vector<std::string>& field_names, 
                                                  const Eigen::VectorXd& data_x,
                                                  const Eigen::VectorXd& data_y) {

    namespace fs = std::filesystem;
    fs::path p(filename);

    EnsureDir ensureDir;
    if (p.has_parent_path()) {
        ensureDir.prepareOutputDirectory(p.parent_path().string());
    }

    std::ofstream out(filename);
    out.imbue(std::locale("C")); 

    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    if (data_x.size() != data_y.size()) {
        Logger::error("Les vecteurs de données X et Y n'ont pas la même taille.");
        return;
    }

    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Vector Field\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    auto nodes = m_mesh.get_nodes();
    out << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << "\n";
    }

    if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        auto elements = m_mesh.get_elements();
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 4) << "\n";
        for (const auto& tri : elements) {
            out << "3 " << tri.node_indices[0] << " " << tri.node_indices[1] << " " << tri.node_indices[2] << "\n";
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "5\n"; 
        }
    } 
    else if (m_mesh.get_element_type() == ElementType::QUAD4) {
        auto elements = m_mesh.get_elements();
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 5) << "\n";
        for (const auto& quad : elements) {
            out << "4 " << quad.node_indices[0] << " " << quad.node_indices[1] << " " << quad.node_indices[2] << " " << quad.node_indices[3] << "\n";
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "9\n"; 
        }
    } 
    else {
        Logger::error("Type d'élément non pris en charge pour l'export VTK.");
        return;
    }

    if (data_x.size() != nodes.size()) {
        Logger::error("Incohérence VTK : Le nombre de vecteurs (" + std::to_string(data_x.size()) + 
                      ") ne correspond pas au nombre de noeuds (" + std::to_string(nodes.size()) + ").");
        return;
    }

    out << "\nPOINT_DATA " << nodes.size() << "\n";
    std::string vector_name = field_names.empty() ? "VectorField" : field_names[0];
    out << "VECTORS " << vector_name << " double\n";
    
    for (size_t j = 0; j < data_x.size(); ++j) {
        out << data_x[j] << " " << data_y[j] << " 0.0\n";
    }

    out.close();
}

inline double swap_double(double val) {
    union {
        double d;
        uint64_t i;
    } u;
    u.d = val;
    u.i = ((u.i & 0x00000000000000FFULL) << 56) |
          ((u.i & 0x000000000000FF00ULL) << 40) |
          ((u.i & 0x0000000000FF0000ULL) << 24) |
          ((u.i & 0x00000000FF000000ULL) << 8)  |
          ((u.i & 0x000000FF00000000ULL) >> 8)  |
          ((u.i & 0x0000FF0000000000ULL) >> 24) |
          ((u.i & 0x00FF000000000000ULL) >> 40) |
          ((u.i & 0xFF00000000000000ULL) >> 56);
    return u.d;
}

void ResultsExporter::exportMultiPhysicsVTK(const std::string& filename, 
                                            const std::vector<std::string>& scalar_names, 
                                            const std::vector<const Eigen::VectorXd*>& scalar_fields,
                                            const std::vector<std::string>& vector_names, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_x, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_y) {

    namespace fs = std::filesystem;
    fs::path p(filename);

    EnsureDir ensureDir;
    if (p.has_parent_path()) {
        ensureDir.prepareOutputDirectory(p.parent_path().string());
    }

    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // --- HELPER POUR LE BINAIRE (Big-Endian) ---
    auto swap_int32 = [](int32_t val) -> int32_t {
        return ((val & 0x000000FF) << 24) |
               ((val & 0x0000FF00) <<  8) |
               ((val & 0x00FF0000) >>  8) |
               ((val & 0xFF000000) >> 24);
    };

    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Multi-Physics\nBINARY\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    // =========================================================
    // 1. GÉOMÉTRIE (POINTS) EN BINAIRE
    // =========================================================
    auto nodes = m_mesh.get_nodes();
    out << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        double x_swap = swap_double(node.x);
        double y_swap = swap_double(node.y);
        out.write(reinterpret_cast<const char*>(&x_swap), sizeof(double));
        out.write(reinterpret_cast<const char*>(&y_swap), sizeof(double));
    }
    out << "\n"; 

    // =========================================================
    // 2. CONNECTIVITÉ (CELLS & CELL_TYPES) EN BINAIRE
    // =========================================================
    
    if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        auto elements = m_mesh.get_elements();
        out << "CELLS " << elements.size() << " " << (elements.size() * 4) << "\n";
        
        for (const auto& tri : elements) {
            int32_t size_swap = swap_int32(3);
            int32_t n1_swap = swap_int32(tri.node_indices[0]);
            int32_t n2_swap = swap_int32(tri.node_indices[1]);
            int32_t n3_swap = swap_int32(tri.node_indices[2]);
            
            out.write(reinterpret_cast<const char*>(&size_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n1_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n2_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n3_swap), sizeof(int32_t));
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            int32_t type_swap = swap_int32(5); 
            out.write(reinterpret_cast<const char*>(&type_swap), sizeof(int32_t));
        }
        out << "\n";
    } 
    else if (m_mesh.get_element_type() == ElementType::QUAD4) {
        auto elements = m_mesh.get_elements();
        out << "CELLS " << elements.size() << " " << (elements.size() * 5) << "\n";
        
        for (const auto& quad : elements) {
            int32_t size_swap = swap_int32(4);
            int32_t n1_swap = swap_int32(quad.node_indices[0]);
            int32_t n2_swap = swap_int32(quad.node_indices[1]);
            int32_t n3_swap = swap_int32(quad.node_indices[2]);
            int32_t n4_swap = swap_int32(quad.node_indices[3]);
            
            out.write(reinterpret_cast<const char*>(&size_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n1_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n2_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n3_swap), sizeof(int32_t));
            out.write(reinterpret_cast<const char*>(&n4_swap), sizeof(int32_t));
        }
        
        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            int32_t type_swap = swap_int32(9); 
            out.write(reinterpret_cast<const char*>(&type_swap), sizeof(int32_t));
        }
        out << "\n";
    }

    // =========================================================
    // 3. DONNÉES PHYSIQUES (POINT_DATA)
    // =========================================================
    out << "POINT_DATA " << nodes.size() << "\n";

    for (size_t i = 0; i < scalar_names.size(); ++i) {
        
        if (scalar_fields[i]->size() != nodes.size()) continue; 

        out << "SCALARS " << scalar_names[i] << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        
        for (const auto& val : *scalar_fields[i]) {
            double swapped = swap_double(val);
            out.write(reinterpret_cast<const char*>(&swapped), sizeof(double));
        }
        out << "\n";
    }

    for (size_t i = 0; i < vector_names.size(); ++i) {
        const Eigen::VectorXd& data_x = *vector_fields_x[i];
        const Eigen::VectorXd& data_y = *vector_fields_y[i];
        
        
        if (data_x.size() != nodes.size() || data_y.size() != nodes.size()) continue; 

        out << "VECTORS " << vector_names[i] << " double\n";
        
        double z_zero = 0.0;
        double z_swapped = swap_double(z_zero);

        for (size_t j = 0; j < data_x.size(); ++j) {
            double x_swapped = swap_double(data_x[j]);
            double y_swapped = swap_double(data_y[j]);
            
            out.write(reinterpret_cast<const char*>(&x_swapped), sizeof(double));
            out.write(reinterpret_cast<const char*>(&y_swapped), sizeof(double));
            out.write(reinterpret_cast<const char*>(&z_swapped), sizeof(double));
        }
        out << "\n";
    }

    out.close();
}

void ResultsExporter::exportConvergenceData(const std::string& filename, 
                                            const std::vector<std::string>& column_names, 
                                            const std::vector<std::vector<double>>& data_columns) {
    // Vérification et création du dossier
    EnsureDir ensureDir;
    ensureDir.prepareOutputDirectory(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // Écriture de l'en-tête
    for (size_t i = 0; i < column_names.size(); ++i) {
        out << column_names[i];
        if (i < column_names.size() - 1) out << "\t";
    }
    out << "\n";

    // Écriture des données ligne par ligne
    size_t num_rows = data_columns[0].size();
    for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < data_columns.size(); ++col) {
            out << data_columns[col][row];
            if (col < data_columns.size() - 1) out << "\t";
        }
        out << "\n";
    }
    out.close();
}