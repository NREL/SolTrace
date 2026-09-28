#pragma once

#include "analysis/rays/ray_geometry.h"
#include "scene_models/element_models.h"
#include "data/simulationresult.h"
#include "support/qt_helpers.h"

#include <QObject>
#include <QVector3D>

namespace SolTrace::GUI::App {

/**
 * @class IntersectionsModule
 * @brief Ray intersection analysis module.
 *
 * Provides access to intersection results from the simulation.
 * Owns the RayGeometry adapter used by the 3D viewport and points it at the
 * currently selected simulation result.
 *
 * QML access pattern: App.intersections.results
 */
class IntersectionsModule : public QObject {
    Q_OBJECT

    SolTrace::GUI::Data::SimulationResultPtr m_results;

    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::AllElementsModel, entity_model)
    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Analysis::RayGeometry, ray_geometry)

public:
    explicit IntersectionsModule(QObject* parent = nullptr);

public slots:
    /// Set the result set whose rays should be visualized.
    void set_results(SolTrace::GUI::Data::SimulationResultPtr);
};
} // namespace SolTrace::GUI::App
