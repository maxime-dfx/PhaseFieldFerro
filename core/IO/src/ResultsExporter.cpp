#include <fstream>
#include <iostream>
#include <cstdint>
#include <filesystem>
#include <vector>
#include "IO/include/ResultsExporter.h"
#include "Utils/include/Profiling.h"
#include "Utils/include/Logger.h"
#include "Utils/include/EnsureDir.h" 
#include <Eigen/Dense>

#ifdef _OPENMP
#include <omp.h>
#endif

ResultsExporter::ResultsExporter(const Mesh& msh) : m_mesh(msh) {}

namespace fs = std::filesystem;

void ResultsExporter::exportToMesh(const std::string& filename) {
    PROFILE_ZONE_NC("ResultsExporter::exportToMesh", PROFILE_COLOR_SEQUENTIAL);
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

    out << "MeshVersionFormatted 1\nDimension 2\n\n";
    
    auto nodes = m_mesh.get_nodes();
    out << "Vertices\n" << nodes.size() << "\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << " "  << node.ref_tag << "\n";
    }
    out << "\n";

    if (m_mesh.get_element_type() == ElementType::QUAD4) {
        auto quads = m_mesh.get_elements();
        out << "Quadrilaterals\n" << quads.size() << "\n";
        for (int i = 0; i < static_cast<int>(quads.size()); ++i) {
            out << m_mesh.get_node_index(i, 0) << " " << m_mesh.get_node_index(i, 1) << " " 
                << m_mesh.get_node_index(i, 2) << " " << m_mesh.get_node_index(i, 3) << " " << quads[i].ref_tag << "\n";
        }
    } else if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        auto tris = m_mesh.get_elements();
        out << "Triangles\n" << tris.size() << "\n";
        for (int i = 0; i < static_cast<int>(tris.size()); ++i) {
            out << (m_mesh.get_node_index(i, 0) + 1) << " " << (m_mesh.get_node_index(i, 1) + 1) << " " 
                << (m_mesh.get_node_index(i, 2) + 1) << " " << tris[i].ref_tag << "\n";        
        }
    } else {
        Logger::error("Type d'élément non pris en charge.");
    }
    
    out << "\nEnd\n";
    out.close();
}

void ResultsExporter::exportCellScalarVTK(const std::string& filename,
                                          const std::string& field_name,
                                          const Eigen::VectorXd& data_field) {
    PROFILE_ZONE_NC("ResultsExporter::exportCellScalarVTK", PROFILE_COLOR_SEQUENTIAL);
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

    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro SimulationManager\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    auto nodes = m_mesh.get_nodes();
    out << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << " 0.0\n";
    }

    auto elements = m_mesh.get_elements();

    if (data_field.size() != static_cast<Eigen::Index>(elements.size())) {
        Logger::error("La taille du champ cellule '", field_name, "' ne correspond pas au nombre d'elements.");
        return;
    }

    if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 4) << "\n";
        for (int i = 0; i < static_cast<int>(elements.size()); ++i) {
            out << "3 " << m_mesh.get_node_index(i, 0) << " " << m_mesh.get_node_index(i, 1) << " " << m_mesh.get_node_index(i, 2) << "\n";
        }

        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "5\n";
        }
    }
    else if (m_mesh.get_element_type() == ElementType::QUAD4) {
        out << "\nCELLS " << elements.size() << " " << (elements.size() * 5) << "\n";
        for (int i = 0; i < static_cast<int>(elements.size()); ++i) {
            out << "4 " << m_mesh.get_node_index(i, 0) << " " << m_mesh.get_node_index(i, 1) << " " << m_mesh.get_node_index(i, 2) << " " << m_mesh.get_node_index(i, 3) << "\n";
        }

        out << "\nCELL_TYPES " << elements.size() << "\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            out << "9\n";
        }
    }
    else {
        Logger::error("Type d'element non pris en charge pour l'export VTK.");
        return;
    }

    out << "\nCELL_DATA " << elements.size() << "\n";
    out << "SCALARS " << field_name << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& val : data_field) {
        out << val << "\n";
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

inline int32_t swap_int32(int32_t val) {
    return ((val & 0x000000FF) << 24) |
           ((val & 0x0000FF00) <<  8) |
           ((val & 0x00FF0000) >>  8) |
           ((val & 0xFF000000) >> 24);
}

void ResultsExporter::exportMultiPhysicsVTK(const std::string& filename, 
                                            const std::vector<std::string>& scalar_names, 
                                            const std::vector<const Eigen::VectorXd*>& scalar_fields,
                                            const std::vector<std::string>& vector_names, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_x, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_y) {
    PROFILE_ZONE_NC("ResultsExporter::exportMultiPhysicsVTK", PROFILE_COLOR_SEQUENTIAL);
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

    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Multi-Physics\nBINARY\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    auto nodes = m_mesh.get_nodes();
    auto elements = m_mesh.get_elements();
    int num_nodes = static_cast<int>(nodes.size());
    int num_elements = static_cast<int>(elements.size());

    // =========================================================
    // 1. POINTS (Préparation OpenMP + Écriture Bloc)
    // =========================================================
    out << "POINTS " << num_nodes << " double\n";
    std::vector<double> points_buffer(num_nodes * 3);
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < num_nodes; ++i) {
        points_buffer[i * 3 + 0] = swap_double(nodes[i].x);
        points_buffer[i * 3 + 1] = swap_double(nodes[i].y);
        points_buffer[i * 3 + 2] = swap_double(0.0);
    }
    out.write(reinterpret_cast<const char*>(points_buffer.data()), points_buffer.size() * sizeof(double));
    out << "\n"; 

    // =========================================================
    // 2. CELLS & CELL_TYPES
    // =========================================================
    if (m_mesh.get_element_type() == ElementType::TRIANGLE3) {
        out << "CELLS " << num_elements << " " << (num_elements * 4) << "\n";
        std::vector<int32_t> cells_buffer(num_elements * 4);
        
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_elements; ++i) {
            cells_buffer[i * 4 + 0] = swap_int32(3);
            cells_buffer[i * 4 + 1] = swap_int32(m_mesh.get_node_index(i, 0));
            cells_buffer[i * 4 + 2] = swap_int32(m_mesh.get_node_index(i, 1));
            cells_buffer[i * 4 + 3] = swap_int32(m_mesh.get_node_index(i, 2));
        }
        out.write(reinterpret_cast<const char*>(cells_buffer.data()), cells_buffer.size() * sizeof(int32_t));
        
        out << "\nCELL_TYPES " << num_elements << "\n";
        std::vector<int32_t> types_buffer(num_elements);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_elements; ++i) {
            types_buffer[i] = swap_int32(5); 
        }
        out.write(reinterpret_cast<const char*>(types_buffer.data()), types_buffer.size() * sizeof(int32_t));
        out << "\n";
    } 
    else if (m_mesh.get_element_type() == ElementType::QUAD4) {
        out << "CELLS " << num_elements << " " << (num_elements * 5) << "\n";
        std::vector<int32_t> cells_buffer(num_elements * 5);
        
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_elements; ++i) {
            cells_buffer[i * 5 + 0] = swap_int32(4);
            cells_buffer[i * 5 + 1] = swap_int32(m_mesh.get_node_index(i, 0));
            cells_buffer[i * 5 + 2] = swap_int32(m_mesh.get_node_index(i, 1));
            cells_buffer[i * 5 + 3] = swap_int32(m_mesh.get_node_index(i, 2));
            cells_buffer[i * 5 + 4] = swap_int32(m_mesh.get_node_index(i, 3));
        }
        out.write(reinterpret_cast<const char*>(cells_buffer.data()), cells_buffer.size() * sizeof(int32_t));
        
        out << "\nCELL_TYPES " << num_elements << "\n";
        std::vector<int32_t> types_buffer(num_elements);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_elements; ++i) {
            types_buffer[i] = swap_int32(9); 
        }
        out.write(reinterpret_cast<const char*>(types_buffer.data()), types_buffer.size() * sizeof(int32_t));
        out << "\n";
    }

    // =========================================================
    // 3. DONNÉES PHYSIQUES POINT_DATA
    // =========================================================
    out << "POINT_DATA " << num_nodes << "\n";

    // EXPORT: Champ d'initialisation des NOEUDS (ref_tag)
    out << "SCALARS node_ref_tag int 1\nLOOKUP_TABLE default\n";
    std::vector<int32_t> node_tag_buffer(num_nodes);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < num_nodes; ++i) {
        node_tag_buffer[i] = swap_int32(nodes[i].ref_tag);
    }
    out.write(reinterpret_cast<const char*>(node_tag_buffer.data()), node_tag_buffer.size() * sizeof(int32_t));
    out << "\n";

    // Buffer réutilisable pour les scalaires
    std::vector<double> point_scalar_buffer(num_nodes);
    
    for (size_t i = 0; i < scalar_names.size(); ++i) {
        if (static_cast<int>(scalar_fields[i]->size()) != num_nodes) continue; 

        out << "SCALARS " << scalar_names[i] << " double 1\nLOOKUP_TABLE default\n";
        const Eigen::VectorXd& current_field = *scalar_fields[i];
        
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < num_nodes; ++j) {
            point_scalar_buffer[j] = swap_double(current_field[j]);
        }
        out.write(reinterpret_cast<const char*>(point_scalar_buffer.data()), point_scalar_buffer.size() * sizeof(double));
        out << "\n";
    }

    // Buffer réutilisable pour les vecteurs
    std::vector<double> point_vector_buffer(num_nodes * 3);
    for (size_t i = 0; i < vector_names.size(); ++i) {
        const Eigen::VectorXd& data_x = *vector_fields_x[i];
        const Eigen::VectorXd& data_y = *vector_fields_y[i];
        
        if (static_cast<int>(data_x.size()) != num_nodes || static_cast<int>(data_y.size()) != num_nodes) continue; 

        out << "VECTORS " << vector_names[i] << " double\n";
        
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < num_nodes; ++j) {
            point_vector_buffer[j * 3 + 0] = swap_double(data_x[j]);
            point_vector_buffer[j * 3 + 1] = swap_double(data_y[j]);
            point_vector_buffer[j * 3 + 2] = swap_double(0.0);
        }
        out.write(reinterpret_cast<const char*>(point_vector_buffer.data()), point_vector_buffer.size() * sizeof(double));
        out << "\n";
    }

    // =========================================================
    // 4. DONNÉES PHYSIQUES CELL_DATA
    // =========================================================
    std::vector<size_t> cell_scalars, cell_vectors;
    for (size_t i = 0; i < scalar_names.size(); ++i) {
        if (static_cast<int>(scalar_fields[i]->size()) == num_elements) cell_scalars.push_back(i);
    }
    for (size_t i = 0; i < vector_names.size(); ++i) {
        if (static_cast<int>(vector_fields_x[i]->size()) == num_elements && static_cast<int>(vector_fields_y[i]->size()) == num_elements) cell_vectors.push_back(i);
    }

    out << "CELL_DATA " << num_elements << "\n";

    // EXPORT: Champ d'initialisation des CELLULES (ref_tag)
    out << "SCALARS cell_ref_tag int 1\nLOOKUP_TABLE default\n";
    std::vector<int32_t> cell_tag_buffer(num_elements);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < num_elements; ++i) {
        cell_tag_buffer[i] = swap_int32(elements[i].ref_tag);
    }
    out.write(reinterpret_cast<const char*>(cell_tag_buffer.data()), cell_tag_buffer.size() * sizeof(int32_t));
    out << "\n";

    if (!cell_scalars.empty() || !cell_vectors.empty()) {
        std::vector<double> cell_scalar_buffer(num_elements);
        for (size_t idx : cell_scalars) {
            out << "SCALARS " << scalar_names[idx] << " double 1\nLOOKUP_TABLE default\n";
            const Eigen::VectorXd& current_field = *scalar_fields[idx];
            
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < num_elements; ++j) {
                cell_scalar_buffer[j] = swap_double(current_field[j]);
            }
            out.write(reinterpret_cast<const char*>(cell_scalar_buffer.data()), cell_scalar_buffer.size() * sizeof(double));
            out << "\n";
        }

        std::vector<double> cell_vector_buffer(num_elements * 3);
        for (size_t idx : cell_vectors) {
            const Eigen::VectorXd& data_x = *vector_fields_x[idx];
            const Eigen::VectorXd& data_y = *vector_fields_y[idx];
            
            out << "VECTORS " << vector_names[idx] << " double\n";
            
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < num_elements; ++j) {
                cell_vector_buffer[j * 3 + 0] = swap_double(data_x[j]);
                cell_vector_buffer[j * 3 + 1] = swap_double(data_y[j]);
                cell_vector_buffer[j * 3 + 2] = swap_double(0.0);
            }
            out.write(reinterpret_cast<const char*>(cell_vector_buffer.data()), cell_vector_buffer.size() * sizeof(double));
            out << "\n";
        }
    }

    out.close();
}