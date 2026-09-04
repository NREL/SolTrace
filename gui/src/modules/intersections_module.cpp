#include "intersections_module.h"

#include <QDebug>

namespace SolTrace::GUI::App {

IntersectionsModule::IntersectionsModule(QObject* parent)
    : QObject(parent),
      m_entity_model(new db::AllElementsModel(this)),
      m_ray_geometry(new analysis::RayGeometry) {
    m_ray_geometry->setParent(this);
}

void IntersectionsModule::set_results(db::SimulationResultPtr ptr) {
    m_results = ptr;
    m_entity_model->reset(nullptr);
    m_ray_geometry->set_results(ptr);

    if (!ptr) return;

    m_entity_model->reset(const_cast<db::Database*>(ptr->database.get()));
}


} // namespace SolTrace::GUI::App
