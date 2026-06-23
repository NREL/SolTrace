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

public:
    explicit ScriptDBInterface(db::Database*, QObject* parent = nullptr);

public slots:

    QJsonArray vec3(double);
    QJsonArray vec3(double, double, double);
    QJsonArray vec3_add(QJsonValue, QJsonValue);
    QJsonArray vec3_sub(QJsonValue, QJsonValue);
    QJsonArray vec3_scale(QJsonValue, double);
    double     vec3_dot(QJsonValue, QJsonValue);
    QJsonArray vec3_cross(QJsonValue, QJsonValue);
    double     vec3_length(QJsonValue);
    double     vec3_distance(QJsonValue, QJsonValue);
    QJsonArray vec3_normalize(QJsonValue);

    QJsonArray quat(double, double, double, double);
    QJsonArray quat_identity();
    QJsonArray quat_from_axis_angle(QJsonValue, double);
    QJsonArray quat_mul(QJsonValue, QJsonValue);
    QJsonArray quat_conjugate(QJsonValue);
    QJsonArray quat_inverse(QJsonValue);
    QJsonArray quat_normalize(QJsonValue);
    QJsonArray quat_rotate_vec3(QJsonValue, QJsonValue);

    db::Entity create();
    void       destroy(db::Entity);
    bool       valid(db::Entity);

    QString get_identity(db::Entity);
    void    set_identity(db::Entity, QString);

    bool get_invisible(db::Entity);
    void set_invisible(db::Entity, bool);

    QJsonObject get_transform(db::Entity);
    void        set_transform(db::Entity, QJsonObject);

    QVector<db::Entity> get_all_materials();
    db::Entity          create_material();
    QJsonObject         get_material_properties(db::Entity);
    void                set_material_properties(db::Entity, QJsonObject);
    void                remove_material(db::Entity);

    QVector<db::Entity> get_all_geometries();
    db::Entity          create_geometry();
    QJsonObject         get_geometry_properties(db::Entity);
    void                set_geometry_properties(db::Entity, QJsonObject);
    void                remove_geometry(db::Entity);

    db::Entity get_material_of(db::Entity);
    void       set_material_of(db::Entity, db::Entity);

    db::Entity get_geometry_of(db::Entity);
    void       set_geometry_of(db::Entity, db::Entity);

signals:
};

} // namespace SolTrace::GUI::Script
