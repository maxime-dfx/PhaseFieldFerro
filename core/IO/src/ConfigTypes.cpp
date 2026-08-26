#include "IO/include/ConfigTypes.h"
#include <stdexcept>

MaterialType material_type_from_string(const std::string& s) {
    if (s == "Ferroelectric")  return MaterialType::Ferroelectric;
    if (s == "PureElastic")    return MaterialType::PureElastic;
    if (s == "PolymerElastic") return MaterialType::PolymerElastic;
    if (s == "GrainBoundary")  return MaterialType::GrainBoundary;
    throw std::invalid_argument(
        "type de matériau inconnu : \"" + s + "\" (attendu : Ferroelectric, "
        "PureElastic, PolymerElastic, GrainBoundary)");
}

std::string material_type_to_string(MaterialType type) {
    switch (type) {
        case MaterialType::Ferroelectric:  return "Ferroelectric";
        case MaterialType::PureElastic:    return "PureElastic";
        case MaterialType::PolymerElastic: return "PolymerElastic";
        case MaterialType::GrainBoundary:  return "GrainBoundary";
    }
    throw std::invalid_argument("MaterialType invalide (valeur enum inconnue)");
}
