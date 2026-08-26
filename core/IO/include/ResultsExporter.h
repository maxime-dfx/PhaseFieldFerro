#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Mesh/include/Mesh.h"

// `using namespace std;` retire (audit §8) : polluait tout traducteur
// incluant ce header. Toutes les utilisations ci-dessous sont deja
// qualifiees std:: (aucun changement de code necessaire).

class ResultsExporter
{
    private:
        const Mesh& m_mesh;
        
        void ensure_directory_exists(const std::string& path);

    public:
        ResultsExporter(const Mesh& msh);

        void exportToMesh(const std::string& filename);     

        // Comme exportScalarVTK (desormais retire, non utilise), mais pour un champ CONSTANT PAR ELEMENT
        // (CELL_DATA), pas par noeud (POINT_DATA). Utilise pour les champs
        // qui n'ont physiquement de sens qu'a l'echelle de l'element, comme
        // grain_id (un noeud a la frontiere de deux grains n'a pas un
        // grain_id unique, contrairement a un champ physique continu).
        void exportCellScalarVTK(const std::string& filename,
                                    const std::string& field_name,
                                    const Eigen::VectorXd& data_field);

        void exportMultiPhysicsVTK(const std::string& filename, 
                                            const std::vector<std::string>& scalar_names, 
                                            const std::vector<const Eigen::VectorXd*>& scalar_fields,
                                            const std::vector<std::string>& vector_names, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_x, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_y); 

};