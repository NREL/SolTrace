#pragma once

#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <QObject>
#include <QPointer>
#include <QVector3D>
#include <QVector>

#include "database/database.h"

namespace SolTrace::GUI::Script {

/// JavaScript-facing database API exposed to user scripts as `db`.
///
/// Methods intentionally use JSON-compatible values at the boundary so scripts
/// do not need to construct C++/Qt value types directly.
class ScriptDBInterface : public QObject {
    Q_OBJECT

    QPointer<db::Database> m_database;
    QString                m_working_directory;

public:
    explicit ScriptDBInterface(db::Database*, QObject* parent = nullptr);

    /// Update the base directory used by get_text_content/get_json_content.
    void update_working_directory(QString);

public slots:
    /// Construct a 3-vector with all components set to value.
    QJsonArray vec3(double value);

    /// Construct a 3-vector from explicit components.
    QJsonArray vec3(double x, double y, double z);

    /// Vector addition.
    QJsonArray vec3_add(QJsonValue a, QJsonValue b);

    /// Vector subtraction.
    QJsonArray vec3_sub(QJsonValue a, QJsonValue b);

    /// Scale a vector by a scalar.
    QJsonArray vec3_scale(QJsonValue value, double scale);

    /// Dot product of two vectors.
    double vec3_dot(QJsonValue a, QJsonValue b);

    /// Cross product of two vectors.
    QJsonArray vec3_cross(QJsonValue a, QJsonValue b);

    /// Euclidean length of a vector.
    double vec3_length(QJsonValue value);

    /// Euclidean distance between two vectors.
    double vec3_distance(QJsonValue a, QJsonValue b);

    /// Normalize a vector.
    QJsonArray vec3_normalize(QJsonValue value);

    /// Construct a quaternion in w, x, y, z order.
    QJsonArray quat(double w, double x, double y, double z);

    /// Identity quaternion.
    QJsonArray quat_identity();

    /// Construct a quaternion rotation from axis and degrees.
    QJsonArray quat_from_axis_angle(QJsonValue axis, double degrees);

    /// Quaternion multiplication.
    QJsonArray quat_mul(QJsonValue a, QJsonValue b);

    /// Quaternion conjugate.
    QJsonArray quat_conjugate(QJsonValue value);

    /// Quaternion inverse.
    QJsonArray quat_inverse(QJsonValue value);

    /// Normalize a quaternion.
    QJsonArray quat_normalize(QJsonValue value);

    /// Rotate a vector by a quaternion.
    QJsonArray quat_rotate_vec3(QJsonValue rotation, QJsonValue value);

    /// Read the database-wide ray source as a JSON object.
    QJsonObject get_ray_source();

    /// Replace the database-wide ray source from a JSON object.
    void set_ray_source(QJsonObject source);

    /// Set a directional sun vector.
    void set_sun_direction(QJsonValue direction);

    /// Set a point-source sun position.
    void set_sun_position(QJsonValue position);

    /// Set the sun shape definition.
    void set_sun_shape(QJsonObject shape);

    /// Return all element entities in the database.
    QVector<db::Entity> get_all_elements();

    /// Create a new element entity.
    db::Entity create();

    /// Destroy an entity and its database-owned relationships.
    void destroy(db::Entity entity);

    /// Check whether an entity is valid in the current database.
    bool valid(db::Entity entity);

    /// Get an entity display name.
    QString get_identity(db::Entity entity);

    /// Set an entity display name.
    void set_identity(db::Entity entity, QString name);

    /// Get whether an entity is hidden in the scene.
    bool get_invisible(db::Entity entity);

    /// Set whether an entity is hidden in the scene.
    void set_invisible(db::Entity entity, bool invisible);

    /// Get local transform as JSON.
    QJsonObject get_transform(db::Entity entity);

    /// Set local transform from JSON.
    void set_transform(db::Entity entity, QJsonObject transform);

    /// Return all material group entities.
    QVector<db::Entity> get_all_materials();

    /// Create a material group.
    db::Entity create_material();

    /// Get material optical properties as JSON.
    QJsonObject get_material_properties(db::Entity material);

    /// Set material optical properties from JSON.
    void set_material_properties(db::Entity material, QJsonObject properties);

    /// Remove a material group.
    void remove_material(db::Entity material);

    /// Return all geometry group entities.
    QVector<db::Entity> get_all_geometries();

    /// Create a geometry group.
    db::Entity create_geometry();

    /// Get geometry surface/aperture properties as JSON.
    QJsonObject get_geometry_properties(db::Entity geometry);

    /// Set geometry surface/aperture properties from JSON.
    void set_geometry_properties(db::Entity geometry, QJsonObject properties);

    /// Remove a geometry group.
    void remove_geometry(db::Entity geometry);

    /// Get the material group assigned to an element.
    db::Entity get_material_of(db::Entity entity);

    /// Assign a material group to an element.
    void set_material_of(db::Entity entity, db::Entity material);

    /// Get the geometry group assigned to an element.
    db::Entity get_geometry_of(db::Entity entity);

    /// Assign a geometry group to an element.
    void set_geometry_of(db::Entity entity, db::Entity geometry);

    /// Read a text file relative to the current working directory.
    QString get_text_content(QString relative_path);

    /// Read and parse a JSON file relative to the current working directory.
    QJsonValue get_json_content(QString relative_path);
};

} // namespace SolTrace::GUI::Script
