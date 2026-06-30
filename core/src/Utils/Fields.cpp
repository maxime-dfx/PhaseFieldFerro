// #include "Utils/Fields.h"
// #include "Physics/Mechanics.h"
// #include "Physics/Electrostatics.h"
// #include "Physics/Polarization.h"
// #include "Physics/Fracture.h"
// #include <Eigen/Dense>

// Fields::Fields(const Mesh& m, const Math& ma, const Datafile& c, 
//                const Electrostatics& e, const Polarization& p, 
//                const Mechanics& me, const Fracture& f)
//     : mesh(m), math(ma), config(c), electrostatics(e), 
//       polarization(p), mechanics(me), fracture(f) 
// {
// }

// const Eigen::VectorXd& Fields::get_Px() const {
//     return polarization.get_Px();
// }
// const Eigen::VectorXd& Fields::get_Py() const {
//     return polarization.get_Py();
// }
// const Eigen::VectorXd& Fields::get_Ex() const {
//     return electrostatics.get_Ex();
// }
// const Eigen::VectorXd& Fields::get_Ey() const {
//     return electrostatics.get_Ey();
// }
// const Eigen::VectorXd& Fields::get_v() const {
//     return fracture.get_v();
// }
// const Eigen::VectorXd& Fields::get_ux() const {
//     return mechanics.get_ux();
// }
// const Eigen::VectorXd& Fields::get_uy() const {
//     return mechanics.get_uy();
// }