#include "IO/Datafile.h"
#include "Utils/Logger.h" 

// Datafile::Datafile(const std::string& filename) {
//     try {
//         m_config = toml::parse_file(filename);
//     } catch (const toml::parse_error& err) {
//         Logger::error("Erreur lors de la lecture du fichier TOML : " + std::string(err.what()));
//         exit(1);
//     }
// }

std::vector<BoundaryRuleConfig> Datafile::get_boundary_rules() const {
    std::vector<BoundaryRuleConfig> rules;
    
    if (auto arr = m_config["boundary_rule"].as_array()) {
        for (auto& node : *arr) {
            if (auto tbl = node.as_table()) {
                BoundaryRuleConfig r;
                r.field     = (*tbl)["field"].value_or<std::string>("");
                r.bc_type   = (*tbl)["bc_type"].value_or<std::string>("DIRICHLET");
                
                r.shape     = (*tbl)["shape"].value_or<std::string>("edge");
                r.edge_name = (*tbl)["edge_name"].value_or<std::string>("");
                r.cx        = (*tbl)["cx"].value_or<double>(0.0);
                r.cy        = (*tbl)["cy"].value_or<double>(0.0);
                r.radius    = (*tbl)["radius"].value_or<double>(0.0);
                r.px        = (*tbl)["px"].value_or<double>(0.0);
                r.py        = (*tbl)["py"].value_or<double>(0.0);
                r.xmin      = (*tbl)["xmin"].value_or<double>(0.0);
                r.xmax      = (*tbl)["xmax"].value_or<double>(0.0);
                r.ymin      = (*tbl)["ymin"].value_or<double>(0.0);
                r.ymax      = (*tbl)["ymax"].value_or<double>(0.0);

                r.profile   = (*tbl)["profile"].value_or<std::string>("constant");
                r.val       = (*tbl)["val"].value_or<double>(0.0);
                r.val_start = (*tbl)["val_start"].value_or<double>(0.0);
                r.val_end   = (*tbl)["val_end"].value_or<double>(0.0);
                r.t_end     = (*tbl)["t_end"].value_or<double>(1.0);
                r.amplitude = (*tbl)["amplitude"].value_or<double>(0.0);
                r.frequency = (*tbl)["frequency"].value_or<double>(0.0);
                r.offset    = (*tbl)["offset"].value_or<double>(0.0);
                r.P0        = (*tbl)["P0"].value_or<double>(0.0);
                r.x0        = (*tbl)["x0"].value_or<double>(0.0);
                r.epsilon   = (*tbl)["epsilon"].value_or<double>(1e-3);
                r.max_val   = (*tbl)["max_val"].value_or<double>(0.0);
                r.y_center  = (*tbl)["y_center"].value_or<double>(0.0);
                r.width     = (*tbl)["width"].value_or<double>(0.0);
                r.k         = (*tbl)["k"].value_or<double>(0.0);
                r.omega     = (*tbl)["omega"].value_or<double>(0.0);
                
                rules.push_back(r);
            }
        }
    }
    return rules;
}