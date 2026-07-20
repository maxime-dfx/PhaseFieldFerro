#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Mesh/Mesh.h"

using namespace std;

class ResultsExporter
{
    private:
        const Mesh& m_mesh;
        
        void ensure_directory_exists(const std::string& path);

    public:
        ResultsExporter(const Mesh& msh);

        void exportToMesh(const std::string& filename);     

        void exportScalarVTK(const std::string& filename, 
                                    const std::string& field_name, 
                                    const Eigen::VectorXd& data_field);

        void exportVectorVTK(const std::string& filename,
                                    const std::vector<std::string>& field_names,
                                    const Eigen::VectorXd& data_x,
                                    const Eigen::VectorXd& data_y);

        void exportMultiPhysicsVTK(const std::string& filename, 
                                            const std::vector<std::string>& scalar_names, 
                                            const std::vector<const Eigen::VectorXd*>& scalar_fields,
                                            const std::vector<std::string>& vector_names, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_x, 
                                            const std::vector<const Eigen::VectorXd*>& vector_fields_y); 

        void exportConvergenceData(const std::string& filename, 
                                    const std::vector<std::string>& column_names, 
                                    const std::vector<std::vector<double>>& data_columns);

};