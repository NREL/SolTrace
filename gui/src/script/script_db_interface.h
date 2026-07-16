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

class ScriptDBInterface : public QObject {
    Q_OBJECT

    QPointer<db::Database> m_database;
    QString                m_working_directory;

public:
    explicit ScriptDBInterface(db::Database*, QObject* parent = nullptr);

    void update_working_directory(QString);

public slots:
    QJsonArray vec3(double value);
    QJsonArray vec3(double x, double y, double z);
    QJsonArray vec3_add(QJsonValue a, QJsonValue b);
    QJsonArray vec3_sub(QJsonValue a, QJsonValue b);
    QJsonArray vec3_scale(QJsonValue value, double scale);
    double     vec3_dot(QJsonValue a, QJsonValue b);
    QJsonArray vec3_cross(QJsonValue a, QJsonValue b);
    double     vec3_length(QJsonValue value);
    double     vec3_distance(QJsonValue a, QJsonValue b);
    QJsonArray vec3_normalize(QJsonValue value);

    QJsonArray quat(double w, double x, double y, double z);
    QJsonArray quat_identity();
    QJsonArray quat_from_axis_angle(QJsonValue axis, double degrees);
    QJsonArray quat_mul(QJsonValue a, QJsonValue b);
    QJsonArray quat_conjugate(QJsonValue value);
    QJsonArray quat_inverse(QJsonValue value);
    QJsonArray quat_normalize(QJsonValue value);
    QJsonArray quat_rotate_vec3(QJsonValue rotation, QJsonValue value);

    QJsonObject get_ray_source();
    void        set_ray_source(QJsonObject source);
    void        set_sun_direction(QJsonValue direction);
    void        set_sun_position(QJsonValue position);
    void        set_sun_shape(QJsonObject shape);

    QVector<db::Entity> get_all_elements();
    db::Entity create();
    void       destroy(db::Entity entity);
    bool       valid(db::Entity entity);

    QString get_identity(db::Entity entity);
    void    set_identity(db::Entity entity, QString name);

    bool get_invisible(db::Entity entity);
    void set_invisible(db::Entity entity, bool invisible);

    QJsonObject get_transform(db::Entity entity);
    void        set_transform(db::Entity entity, QJsonObject transform);

    QVector<db::Entity> get_all_materials();
    db::Entity          create_material();
    QJsonObject         get_material_properties(db::Entity material);
    void                set_material_properties(db::Entity  material,
                                                QJsonObject properties);
    void                remove_material(db::Entity material);

    QVector<db::Entity> get_all_geometries();
    db::Entity          create_geometry();
    QJsonObject         get_geometry_properties(db::Entity geometry);
    void                set_geometry_properties(db::Entity  geometry,
                                                QJsonObject properties);
    void                remove_geometry(db::Entity geometry);

    db::Entity get_material_of(db::Entity entity);
    void       set_material_of(db::Entity entity, db::Entity material);

    db::Entity get_geometry_of(db::Entity entity);
    void       set_geometry_of(db::Entity entity, db::Entity geometry);

    QString    get_text_content(QString relative_path);
    QJsonValue get_json_content(QString relative_path);

signals:
};

} // namespace SolTrace::GUI::Script
