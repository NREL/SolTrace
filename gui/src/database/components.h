#pragma once

#include <QColor>
#include <QImage>
#include <QObject>
#include <QVector>

#include <entt/entt.hpp>
#include <glm/gtc/quaternion.hpp>

#include "analysis/baked_flux_map.h"
#include "entity.h"
#include "simulation_data_api.hpp"

namespace SD = SolTrace::Data;

namespace db {

// When adding components, be sure to add them to the database clone function!

/// Tag component for visually hiding an element
struct InvisibleComponent { };

/// Tag for disabling an element in the simulation
struct DisabledComponent { };

/// An identity, or non-unique name for an entity (various semantics)
struct IdentityComponent {
    QString name;
};

struct ElementComponent { };

/// Describes the parent of this entity, if it has one. DO NOT modify this
/// component directly!
struct ChildOfComponent {
    Entity parent;

    template <class M>
    void remap_entities(M& mapper) {
        parent = mapper(parent);
    }
};

/// Describes the children of an entity. DO NOT modify this component directly!
struct ChildrenComponent {
    QVector<Entity> children;

    template <class M>
    void remap_entities(M& mapper) {
        for (auto& e : children) {
            e = mapper(e);
        }
    }
};

/// Describe the attitude of this entity.
struct TransformComponent {
    glm::dvec3 position;
    glm::dquat rotation;

    glm::dmat4 as_matrix() const;

    static TransformComponent from_json(QJsonObject const&);
    QJsonObject               to_json() const;
};

/// Describe the global attitude of this entity.
struct GlobalTransformComponent {
    glm::dvec3 position;
    glm::dquat rotation;

    glm::dmat4 as_matrix() const;

    static GlobalTransformComponent compute_for(entt::registry const& reg,
                                                entt::entity          entity);
};

/// A Global describing the ray source.
struct RaySourceResource {
    SD::ray_source_ptr source;

    RaySourceResource clone() const;
};

struct DatabaseNameResource {
    QString name;

    DatabaseNameResource clone() const;
};

/// A set of material properties.
struct MaterialComponent {
    SD::OpticalPropertySet optics;

    bool operator==(MaterialComponent const&) const;
};

/// A group using the same material. DO NOT modify the member information
/// directly.
struct MaterialGroupComponent {
    QVector<Entity> members;

    template <class M>
    void remap_entities(M& mapper) {
        for (auto& e : members) {
            e = mapper(e);
        }
    }
};

/// Describes the material group this entity belongs to. DO NOT modify this
/// component directly!
struct MaterialGroupMemberComponent {
    Entity group;

    template <class M>
    void remap_entities(M& mapper) {
        group = mapper(group);
    }
};

/// A set of geometry properties.
struct GeometryComponent {
    SD::aperture_ptr aperture;
    SD::surface_ptr  surface;

    bool operator==(GeometryComponent const&) const;

    GeometryComponent clone() const;
};

/// A group using the same geometry. DO NOT modify the member information
/// directly.
struct GeometryGroupComponent {
    QVector<Entity> members;

    template <class M>
    void remap_entities(M& mapper) {
        for (auto& e : members) {
            e = mapper(e);
        }
    }
};

/// Describes the geometry group this entity belongs to. DO NOT modify this
/// component directly!
struct GeometryGroupMemberComponent {
    Entity group;

    template <class M>
    void remap_entities(M& mapper) {
        group = mapper(group);
    }
};

/// A component indicating this entity is a Tag description
struct TagComponent { };

/// This is a tag component that is used, with a string, to indicate membership
/// in a tag. Do not modify this directly!
struct ATagMemberComponent { };

/// Lists the string tags this entity has
struct TagMembershipComponent {
    QVector<Entity> tags;

    template <class M>
    void remap_entities(M& mapper) {
        for (auto& e : tags) {
            e = mapper(e);
        }
    }
};

/// Denotes an entity that could not be imported properly
struct ImportErrorComponent {
    QString reason;
};

/// Denotes whether an instance or group has been selected through object
/// picking
struct SelectedComponent { };

/// Color of the model
struct ColorComponent {
    QColor color = Qt::white;
};

/// Marks if this entity has a flux map
/// Should be hidden during instance rendering
struct HasFluxMapComponent {
    analysis::BakedFluxMapPtr map_info;
};

} // namespace db
