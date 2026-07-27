#pragma once

#include "analysis/ray_geometry.h"
#include "database/simulationresult.h"
#include "module_common.h"
#include "utilities/qt_helpers.h"

#include <QObject>
#include <QVector3D>

namespace SolTrace::GUI::App {

/**
 * @class IntersectionsModule
 * @brief Ray intersection analysis module.
 *
 * Provides access to intersection results from the simulation.
 * Holds non-owning references to both the shared results backend
 * and the intersections-specific backend.
 *
 * QML access pattern: App.intersections.results
 */
class IntersectionsModule : public QObject {
    Q_OBJECT

    db::SimulationResultPtr m_results;

    QOBJECT_READONLY_PROPERTY(analysis::RayGeometry, ray_geometry)

public:
    explicit IntersectionsModule(QObject* parent = nullptr);


public slots:
    void set_results(db::SimulationResultPtr);
};


} // namespace SolTrace::GUI::App
