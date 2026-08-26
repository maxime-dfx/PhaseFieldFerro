#pragma once
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "Materials/Core/MaterialConcepts.h"
#include "IO/include/ConfigTypes.h"
#include "Materials/Models/Ferroelectric.h"
#include "Materials/Models/PureElastic.h"
#include "Materials/Models/PolymerElastic.h"
#include "Materials/Models/GrainBoundary.h"

class MaterialManager final {
private:
    std::unordered_map<int, std::unique_ptr<MaterialModel>> m_materials;

public:
    std::unique_ptr<MaterialModel> create(const MaterialConfig& config) {
    switch (config.type) {
        case MaterialType::Ferroelectric:
            return std::make_unique<SingleCrystalMaterial>(config);
        case MaterialType::PureElastic:
            return std::make_unique<PureElasticMaterial>(config);
        case MaterialType::PolymerElastic:
            return std::make_unique<PolymerElasticMaterial>(config);
        case MaterialType::GrainBoundary:
            return std::make_unique<GrainBoundaryMaterial>(config);
    }
    throw std::invalid_argument(
        "MaterialManager::create : MaterialType inconnu pour l'entrée [[materials]] id=" +
        std::to_string(config.id));
    }

    const MaterialModel& create_and_register(int material_id, const MaterialConfig& config) {
        auto material = create(config);
        return *(m_materials.emplace(material_id, std::move(material)).first->second);
    }

    void register_material(int material_id, std::unique_ptr<MaterialModel> material) {
        m_materials[material_id] = std::move(material);
    }

    const MaterialModel& get_material(int material_id) const {
        auto it = m_materials.find(material_id);
        if (it == m_materials.end()) {
            throw std::runtime_error("Materiau ID " + std::to_string(material_id) + " introuvable.");
        }
        return *(it->second);
    }
};
