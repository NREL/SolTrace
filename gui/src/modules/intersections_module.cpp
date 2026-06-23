#include "intersections_module.h"

#include <QDebug>

namespace SolTrace::GUI::App {

IntersectionsModule::IntersectionsModule(QObject* parent)
    : QObject(parent), m_ray_geometry(new analysis::RayGeometry) {
    m_ray_geometry->setParent(this);
}

void IntersectionsModule::set_results(db::SimulationResultPtr ptr) {
    m_results = ptr;
    m_ray_geometry->set_results(ptr);
}


} // namespace SolTrace::GUI::App
