#pragma once

class IDofMap {
public:
    virtual ~IDofMap() = default;
    
    // Récupère l'indice du DDL pour un noeud local d'un élément
    virtual int dof_for(int element_id, int local_node) const = 0;
    
    // Récupère le noeud géométrique global associé à un DDL
    virtual int node_for(int dof_id) const = 0;
    
    // Nombre total de DDLs
    virtual int num_dofs() const = 0;
};