#include "backend.h"

#include "job_control/job_run_common.h"
#include "utilities/math_utility.h"

#include "surface.hpp"

#include <QApplication>
#include <QFileInfo>
#include <QFuture>
#include <QFutureWatcher>
#include <QStringView>
#include <QTextStream>
#include <QtConcurrentRun>

// TODO COORDINATE SYSTEMS. the sim could be arbitrary


ResultsBackend::ResultsBackend(QObject* parent)
    : QObject(parent), m_ray_geometry(new analysis::RayGeometry) { }

void ResultsBackend::set_results(db::SimulationResultPtr ptr) {
    m_ray_geometry->set_results(ptr);
}

analysis::RayGeometry* ResultsBackend::ray_geometry() const {
    return m_ray_geometry.get();
}

//==============================================================================

Backend::Backend(QObject* parent) : QObject(parent) { }


void Backend::install(db::Database* db) {
    qDebug() << Q_FUNC_INFO << db;

    if (m_current_database) { delete m_current_database; }

    m_current_database = db;

    if (db) { db->setParent(this); }

    // m_job_backend->install(db);

    // m_breadcrumb_model->reset(db);
    // m_child_model->reset(db);
    // m_render_groups_model->reset(db);
    // m_tags_model->reset(db);
    // m_instance_edit_model->reset(db);
    // m_world_geometry_model->reset(db);
}

void Backend::reset() { }

void Backend::start_load_file(QUrl file) { }

SunBackend::SunBackend(QObject *parent) :
    QObject(parent)
{}

TracingBackend::TracingBackend(QObject *parent) :
    QObject(parent)
{}

MaterialsBackend::MaterialsBackend(QObject *parent) :
    QObject(parent)
{}

GeometryBackend::GeometryBackend(QObject *parent) :
    QObject(parent)
{}

IntersectionsBackend::IntersectionsBackend(QObject *parent) :
    QObject(parent)
{}

FluxBackend::FluxBackend(QObject *parent) :
    QObject(parent)
{}


void Backend::file_ready()
{
    // stub
}
