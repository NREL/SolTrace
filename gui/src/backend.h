#pragma once

#include "analysis/ray_geometry.h"
#include "database/database.h"
#include "database/database_models.h"
#include "database/worldgeometrymodel.h"
#include "job_control/job_run.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"

#include <QObject>
#include <QQmlEngine>
#include <QSharedPointer>

//==============================================================================

// Shared backend resource, accessed by multiple domain modules (Intersections, Flux)
class ResultsBackend : public QObject {
    Q_OBJECT

    QPointer<db::SimulationResult> m_results;

    std::unique_ptr<analysis::RayGeometry> m_ray_geometry;

    Q_PROPERTY(analysis::RayGeometry* ray_geometry READ ray_geometry NOTIFY
                   ray_geometry_changed FINAL)

public:
    explicit ResultsBackend(QObject* parent = nullptr);
    virtual ~ResultsBackend() = default;

    analysis::RayGeometry* ray_geometry() const;

public slots:
    void set_results(db::SimulationResultPtr);

signals:
    void ray_geometry_changed();
};

//==============================================================================

class SunBackend : public QObject {
    Q_OBJECT

public:
    explicit SunBackend(QObject* parent = nullptr);

    // bridge here
    // translate SolTrace::App::Sun to SolTrace::Data::Sun
};

class TracingBackend : public QObject {
    Q_OBJECT

public:
    explicit TracingBackend(QObject* parent = nullptr);

    // bridge here
    // translate SolTrace::App::Tracing to SolTrace::Data::SimulationData
};

class MaterialsBackend : public QObject {
    Q_OBJECT

public:
    explicit MaterialsBackend(QObject* parent = nullptr);

    QOBJECT_READONLY_PROPERTY(db::BreadcrumbModel, breadcrumb_model);
    QOBJECT_READONLY_PROPERTY(db::ChildModel, child_model);
    QOBJECT_READONLY_PROPERTY(db::MaterialGroupsModel, render_groups_model);
    QOBJECT_READONLY_PROPERTY(db::TagsModel, tags_model);
};

class GeometryBackend : public QObject {
    Q_OBJECT

public:
    explicit GeometryBackend(QObject* parent = nullptr);

    QOBJECT_READONLY_PROPERTY(db::AnInstanceEditor, instance_edit_model);
    QOBJECT_READONLY_PROPERTY(db::WorldGeometryModel, world_geometry_model);
};

class IntersectionsBackend : public QObject {
    Q_OBJECT

public:
    explicit IntersectionsBackend(QObject* parent = nullptr);

    // install db here
};

class FluxBackend : public QObject {
    Q_OBJECT

public:
    explicit FluxBackend(QObject* parent = nullptr);

    // install db here
};

class Backend : public QObject {
    Q_OBJECT

    QML_ELEMENT
    QML_SINGLETON

    // If data came from a file, this is the path
    Q_WRITABLE_PROPERTY(QString, current_data_path, {});

    // Current content
    QPointer<db::Database> m_current_database;

    // Shared resources, such as results
    QOBJECT_READONLY_PROPERTY(ResultsBackend, results);

    // Clear ownership and access semantics using separate modules
    QOBJECT_READONLY_PROPERTY(SunBackend, sun);
    QOBJECT_READONLY_PROPERTY(TracingBackend, tracing);
    QOBJECT_READONLY_PROPERTY(MaterialsBackend, materials);
    QOBJECT_READONLY_PROPERTY(GeometryBackend, geometry);
    QOBJECT_READONLY_PROPERTY(IntersectionsBackend, intersections);
    QOBJECT_READONLY_PROPERTY(FluxBackend, flux);


    void install(db::Database*);

private slots:
    void file_ready();

public:
    explicit Backend(QObject* parent = nullptr);

public slots:
    void reset();
    void start_load_file(QUrl);

signals:
    void notification(ANotification);
};
