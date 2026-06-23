#include "script_db_interface.h"

#include <QDebug>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonValue>

#include <glm/geometric.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

#include <exception>
#include <utility>

namespace SolTrace::GUI::Script {

namespace {

QJsonArray to_json(glm::dvec3 const& v) {
    return QJsonArray { v.x, v.y, v.z };
}

QJsonArray to_json(glm::dquat const& q) {
    return QJsonArray { q.w, q.x, q.y, q.z };
}

double deg_to_rad(double degrees) {
    return degrees * glm::pi<double>() / 180.0;
}

bool read_vec3(QJsonValue const& value, glm::dvec3& out) {
    if (value.isArray()) {
        auto array = value.toArray();
        if (array.size() != 3) return false;
        out = glm::dvec3 {
            array.at(0).toDouble(),
            array.at(1).toDouble(),
            array.at(2).toDouble(),
        };
        return true;
    }

    if (value.isObject()) {
        auto obj = value.toObject();
        out      = glm::dvec3 {
            obj.value("x").toDouble(),
            obj.value("y").toDouble(),
            obj.value("z").toDouble(),
        };
        return obj.contains("x") && obj.contains("y") && obj.contains("z");
    }

    return false;
}

bool read_vec3_or_scalar(QJsonValue const& value, glm::dvec3& out) {
    if (value.isDouble()) {
        auto scalar = value.toDouble();
        out         = glm::dvec3 { scalar, scalar, scalar };
        return true;
    }

    return read_vec3(value, out);
}

bool read_quat(QJsonValue const& value, glm::dquat& out) {
    if (value.isArray()) {
        auto array = value.toArray();
        if (array.size() != 4) return false;
        out = glm::dquat {
            array.at(0).toDouble(),
            array.at(1).toDouble(),
            array.at(2).toDouble(),
            array.at(3).toDouble(),
        };
        return true;
    }

    if (value.isObject()) {
        auto obj = value.toObject();
        out      = glm::dquat {
            obj.value("w").toDouble(),
            obj.value("x").toDouble(),
            obj.value("y").toDouble(),
            obj.value("z").toDouble(),
        };
        return obj.contains("w") && obj.contains("x") && obj.contains("y") &&
               obj.contains("z");
    }

    return false;
}

QJsonObject to_qjson(nlohmann::ordered_json const& json) {
    auto bytes = QByteArray::fromStdString(json.dump());
    return QJsonDocument::fromJson(bytes).object();
}

nlohmann::ordered_json to_njson(QJsonObject const& object) {
    auto bytes = QJsonDocument(object).toJson(QJsonDocument::Compact);
    return nlohmann::ordered_json::parse(bytes.constData());
}

QJsonObject to_qjson(SD::OpticalPropertySet const& properties,
                     SD::OpticalSide               side) {
    double refraction_index_front;
    double refraction_index_back;
    properties.get_refraction_indices(refraction_index_front,
                                      refraction_index_back);

    return QJsonObject {
        { "my_type",
          QString::fromStdString(
              SD::InteractionTypeMap.at(properties.get_interaction_type())) },
        { "error_distribution_type",
          QString::fromStdString(
              SD::DistributionTypeMap.at(
                  properties.get_error_distribution(side))) },
        { "transmissivity", properties.get_transmissivity(side) },
        { "reflectivity", properties.get_reflectivity(side) },
        { "slope_error", properties.get_slope_error(side) },
        { "specularity_error", properties.get_specularity_error(side) },
        { "refraction_index_front", refraction_index_front },
        { "refraction_index_back", refraction_index_back },
    };
}

void patch_optical_properties(SD::OpticalPropertySet& properties,
                              SD::OpticalSide         side,
                              QJsonObject const&      object) {
    if (object.contains("my_type")) {
        auto key = object.value("my_type").toString().toStdString();
        if (auto value = db::reverse_lookup(SD::InteractionTypeMap, key)) {
            properties.set_interaction_type(*value);
        }
    }

    if (object.contains("error_distribution_type")) {
        auto key =
            object.value("error_distribution_type").toString().toStdString();
        if (auto value = db::reverse_lookup(SD::DistributionTypeMap, key)) {
            auto slope_error = properties.get_slope_error(side);
            auto spec_error  = properties.get_specularity_error(side);
            properties.set_errors(side, *value, slope_error, spec_error);
        }
    }

    if (object.contains("transmissivity")) {
        properties.set_transmissivity(
            side, object.value("transmissivity").toDouble());
    }
    if (object.contains("transmitivity")) {
        properties.set_transmissivity(
            side, object.value("transmitivity").toDouble());
    }
    if (object.contains("reflectivity")) {
        properties.set_reflectivity(
            side, object.value("reflectivity").toDouble());
    }
    if (object.contains("slope_error")) {
        auto distribution = properties.get_error_distribution(side);
        auto spec_error   = properties.get_specularity_error(side);
        properties.set_errors(
            side,
            distribution,
            object.value("slope_error").toDouble(),
            spec_error);
    }
    if (object.contains("specularity_error")) {
        auto distribution = properties.get_error_distribution(side);
        auto slope_error  = properties.get_slope_error(side);
        properties.set_errors(
            side,
            distribution,
            slope_error,
            object.value("specularity_error").toDouble());
    }

    double refraction_index_front;
    double refraction_index_back;
    properties.get_refraction_indices(refraction_index_front,
                                      refraction_index_back);

    bool has_refraction_change = false;
    if (object.contains("refraction_index_front")) {
        refraction_index_front =
            object.value("refraction_index_front").toDouble();
        has_refraction_change = true;
    }
    if (object.contains("refraction_index_back")) {
        refraction_index_back =
            object.value("refraction_index_back").toDouble();
        has_refraction_change = true;
    }

    if (has_refraction_change) {
        properties.set_refraction_indices(refraction_index_front,
                                          refraction_index_back);
    }
}

template <class Component>
QVector<db::Entity> collect_entities(entt::registry const& registry) {
    QVector<db::Entity> ret;
    for (auto entity : registry.view<Component const>()) {
        ret.push_back(db::Entity { entity });
    }
    return ret;
}

} // namespace

ScriptDBInterface::ScriptDBInterface(db::Database* database, QObject* parent)
    : QObject { parent }, m_database { database } { }

QJsonArray ScriptDBInterface::vec3(double value) {
    return to_json(glm::dvec3 { value, value, value });
}

QJsonArray ScriptDBInterface::vec3(double x, double y, double z) {
    return to_json(glm::dvec3 { x, y, z });
}

QJsonArray ScriptDBInterface::vec3_add(QJsonValue a, QJsonValue b) {
    glm::dvec3 lhs;
    glm::dvec3 rhs;
    if (!read_vec3_or_scalar(a, lhs) || !read_vec3_or_scalar(b, rhs)) return {};
    return to_json(lhs + rhs);
}

QJsonArray ScriptDBInterface::vec3_sub(QJsonValue a, QJsonValue b) {
    glm::dvec3 lhs;
    glm::dvec3 rhs;
    if (!read_vec3_or_scalar(a, lhs) || !read_vec3_or_scalar(b, rhs)) return {};
    return to_json(lhs - rhs);
}

QJsonArray ScriptDBInterface::vec3_scale(QJsonValue value, double scale) {
    glm::dvec3 vector;
    if (!read_vec3_or_scalar(value, vector)) return {};
    return to_json(vector * scale);
}

double ScriptDBInterface::vec3_dot(QJsonValue a, QJsonValue b) {
    glm::dvec3 lhs;
    glm::dvec3 rhs;
    if (!read_vec3_or_scalar(a, lhs) || !read_vec3_or_scalar(b, rhs)) return 0.0;
    return glm::dot(lhs, rhs);
}

QJsonArray ScriptDBInterface::vec3_cross(QJsonValue a, QJsonValue b) {
    glm::dvec3 lhs;
    glm::dvec3 rhs;
    if (!read_vec3(a, lhs) || !read_vec3(b, rhs)) return {};
    return to_json(glm::cross(lhs, rhs));
}

double ScriptDBInterface::vec3_length(QJsonValue value) {
    glm::dvec3 vector;
    if (!read_vec3_or_scalar(value, vector)) return 0.0;
    return glm::length(vector);
}

double ScriptDBInterface::vec3_distance(QJsonValue a, QJsonValue b) {
    glm::dvec3 lhs;
    glm::dvec3 rhs;
    if (!read_vec3(a, lhs) || !read_vec3(b, rhs)) return 0.0;
    return glm::distance(lhs, rhs);
}

QJsonArray ScriptDBInterface::vec3_normalize(QJsonValue value) {
    glm::dvec3 vector;
    if (!read_vec3(value, vector)) return {};
    auto length = glm::length(vector);
    if (length <= 0.0) return {};
    return to_json(vector / length);
}

QJsonArray ScriptDBInterface::quat(double w, double x, double y, double z) {
    return to_json(glm::dquat { w, x, y, z });
}

QJsonArray ScriptDBInterface::quat_identity() {
    return to_json(glm::dquat { 1.0, 0.0, 0.0, 0.0 });
}

QJsonArray ScriptDBInterface::quat_from_axis_angle(QJsonValue axis,
                                                   double     degrees) {
    glm::dvec3 axis_vector;
    if (!read_vec3(axis, axis_vector)) return {};
    auto length = glm::length(axis_vector);
    if (length <= 0.0) return {};

    return to_json(glm::angleAxis(deg_to_rad(degrees), axis_vector / length));
}

QJsonArray ScriptDBInterface::quat_mul(QJsonValue a, QJsonValue b) {
    glm::dquat lhs;
    glm::dquat rhs;
    if (!read_quat(a, lhs) || !read_quat(b, rhs)) return {};
    return to_json(lhs * rhs);
}

QJsonArray ScriptDBInterface::quat_conjugate(QJsonValue value) {
    glm::dquat q;
    if (!read_quat(value, q)) return {};
    return to_json(glm::conjugate(q));
}

QJsonArray ScriptDBInterface::quat_inverse(QJsonValue value) {
    glm::dquat q;
    if (!read_quat(value, q)) return {};
    auto length = glm::length(q);
    if (length <= 0.0) return {};
    return to_json(glm::inverse(q));
}

QJsonArray ScriptDBInterface::quat_normalize(QJsonValue value) {
    glm::dquat q;
    if (!read_quat(value, q)) return {};
    auto length = glm::length(q);
    if (length <= 0.0) return {};
    return to_json(glm::normalize(q));
}

QJsonArray ScriptDBInterface::quat_rotate_vec3(QJsonValue rotation,
                                               QJsonValue value) {
    glm::dquat q;
    glm::dvec3 vector;
    if (!read_quat(rotation, q) || !read_vec3(value, vector)) return {};
    auto length = glm::length(q);
    if (length <= 0.0) return {};
    return to_json(glm::normalize(q) * vector);
}

db::Entity ScriptDBInterface::create() {
    if (!m_database) return {};

    auto  entity = m_database->create();
    auto& reg    = m_database->as_registry();

    reg.emplace<db::ElementComponent>(entity);
    reg.emplace<db::TransformComponent>(
        entity,
        db::TransformComponent {
            .position = glm::dvec3 { 0.0 },
            .rotation = glm::dquat { 1.0, 0.0, 0.0, 0.0 },
        });

    return db::Entity { entity };
}

void ScriptDBInterface::destroy(db::Entity entity) {
    if (!m_database || !m_database->valid(entity)) return;

    auto& reg = m_database->as_registry();

    if (reg.all_of<db::MaterialGroupComponent>(entity)) {
        m_database->delete_material_group(entity);
        return;
    }

    if (reg.all_of<db::GeometryGroupComponent>(entity)) {
        m_database->delete_geometry_group(entity);
        return;
    }

    if (reg.all_of<db::TagComponent>(entity)) {
        m_database->delete_tag(entity);
        return;
    }

    if (auto children = m_database->children_of(entity); !children.empty()) {
        QVector<entt::entity> copy;
        copy.reserve(children.size());
        for (auto child : children) {
            copy.push_back(child);
        }
        for (auto child : std::as_const(copy)) {
            m_database->unset_parent(child);
        }
    }

    m_database->unset_parent(entity);
    m_database->remove_material(entity);
    m_database->remove_geometry(entity);

    auto                  tags = m_database->tags_for(entity);
    QVector<entt::entity> tag_copy;
    tag_copy.reserve(tags.size());
    for (auto tag : tags) {
        tag_copy.push_back(tag);
    }
    for (auto tag : std::as_const(tag_copy)) {
        m_database->unassign_tag(entity, tag);
    }

    reg.destroy(entity);
}

bool ScriptDBInterface::valid(db::Entity entity) {
    return m_database && m_database->valid(entity);
}

QString ScriptDBInterface::get_identity(db::Entity entity) {
    if (!m_database || !m_database->valid(entity)) return {};
    return m_database->name_of(entity);
}

void ScriptDBInterface::set_identity(db::Entity entity, QString name) {
    if (!m_database || !m_database->valid(entity)) return;
    m_database->identity.set(entity, db::IdentityComponent { .name = name });
}

bool ScriptDBInterface::get_invisible(db::Entity entity) {
    if (!m_database || !m_database->valid(entity)) return false;
    return m_database->as_registry().all_of<db::InvisibleComponent>(entity);
}

void ScriptDBInterface::set_invisible(db::Entity entity, bool invisible) {
    if (!m_database || !m_database->valid(entity)) return;

    if (invisible) {
        m_database->invisible.set(entity, db::InvisibleComponent {});
    } else if (get_invisible(entity)) {
        m_database->invisible.remove(entity);
    }
}

QJsonObject ScriptDBInterface::get_transform(db::Entity entity) {
    if (!m_database) return {};

    auto* transform = m_database->transform.get(entity);
    if (!transform) return {};

    return QJsonObject {
        { "position", to_json(transform->position) },
        { "rotation", to_json(transform->rotation) },
    };
}

void ScriptDBInterface::set_transform(db::Entity entity, QJsonObject object) {
    if (!m_database || !m_database->valid(entity)) return;

    m_database->transform.patch(entity, [&](db::TransformComponent& transform) {
        glm::dvec3 position;
        if (object.contains("position") &&
            read_vec3(object.value("position"), position)) {
            transform.position = position;
        }

        glm::dquat rotation;
        if (object.contains("rotation") &&
            read_quat(object.value("rotation"), rotation)) {
            transform.rotation = glm::normalize(rotation);
        }
    });
}

QVector<db::Entity> ScriptDBInterface::get_all_materials() {
    if (!m_database) return {};
    return collect_entities<db::MaterialGroupComponent>(
        m_database->as_registry());
}

db::Entity ScriptDBInterface::create_material() {
    if (!m_database) return {};
    return db::Entity { m_database->add_material_group("Material", {}) };
}

QJsonObject ScriptDBInterface::get_material_properties(db::Entity entity) {
    if (!m_database) return {};

    auto* material = m_database->material_parameters.get(entity);
    if (!material) return {};

    return QJsonObject {
        { "front", to_qjson(material->optics, SD::OpticalSide::Front) },
        { "back", to_qjson(material->optics, SD::OpticalSide::Back) },
    };
}

void ScriptDBInterface::set_material_properties(db::Entity  entity,
                                                QJsonObject object) {
    if (!m_database || !m_database->valid(entity)) return;

    m_database->material_parameters.try_patch(
        entity, [&](db::MaterialComponent& material) {
            if (object.value("front").isObject()) {
                patch_optical_properties(material.optics,
                                         SD::OpticalSide::Front,
                                         object.value("front").toObject());
            }

            if (object.value("back").isObject()) {
                patch_optical_properties(material.optics,
                                         SD::OpticalSide::Back,
                                         object.value("back").toObject());
            }
        });
}

void ScriptDBInterface::remove_material(db::Entity entity) {
    if (!m_database) return;
    m_database->delete_material_group(entity);
}

QVector<db::Entity> ScriptDBInterface::get_all_geometries() {
    if (!m_database) return {};
    return collect_entities<db::GeometryGroupComponent>(
        m_database->as_registry());
}

db::Entity ScriptDBInterface::create_geometry() {
    if (!m_database) return {};
    return db::Entity { m_database->add_geometry_group("Geometry", {}) };
}

QJsonObject ScriptDBInterface::get_geometry_properties(db::Entity entity) {
    if (!m_database) return {};

    auto* geometry = m_database->geometry_parameters.get(entity);
    if (!geometry) return {};

    QJsonObject ret;

    if (geometry->aperture) {
        nlohmann::ordered_json json;
        geometry->aperture->write_json(json);
        ret["aperture"] = to_qjson(json);
    }

    if (geometry->surface) {
        nlohmann::ordered_json json;
        geometry->surface->write_json(json);
        ret["surface"] = to_qjson(json);
    }

    return ret;
}

void ScriptDBInterface::set_geometry_properties(db::Entity  entity,
                                                QJsonObject object) {
    if (!m_database || !m_database->valid(entity)) return;

    m_database->geometry_parameters.try_patch(
        entity, [&](db::GeometryComponent& geometry) {
            try {
                if (object.value("aperture").isObject()) {
                    geometry.aperture = SD::Aperture::make_aperture_from_json(
                        to_njson(object.value("aperture").toObject()));
                }

                if (object.value("surface").isObject()) {
                    geometry.surface = SD::make_surface_from_json(
                        to_njson(object.value("surface").toObject()));
                }
            } catch (std::exception const& e) {
                qWarning() << "Failed to set geometry properties:" << e.what();
            }
        });
}

void ScriptDBInterface::remove_geometry(db::Entity entity) {
    if (!m_database) return;
    m_database->delete_geometry_group(entity);
}

db::Entity ScriptDBInterface::get_material_of(db::Entity entity) {
    if (!m_database || !m_database->valid(entity)) return {};
    return db::Entity { m_database->material_of(entity) };
}

void ScriptDBInterface::set_material_of(db::Entity entity,
                                        db::Entity material) {
    if (!m_database || !m_database->valid(entity)) return;
    m_database->assign_material(entity, material);
}

db::Entity ScriptDBInterface::get_geometry_of(db::Entity entity) {
    if (!m_database || !m_database->valid(entity)) return {};
    return db::Entity { m_database->geometry_of(entity) };
}

void ScriptDBInterface::set_geometry_of(db::Entity entity,
                                        db::Entity geometry) {
    if (!m_database || !m_database->valid(entity)) return;
    m_database->assign_geometry(entity, geometry);
}

} // namespace SolTrace::GUI::Script
