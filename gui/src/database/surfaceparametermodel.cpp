#include "surfaceparametermodel.h"

#include <magic_enum/magic_enum.hpp>

#include <QtMath>

namespace db {

namespace {

QVector<SurfaceParameter> make_parameters_for(SD::SurfaceType type) {
    QVector<SurfaceParameter> ret;

    switch (type) {
    case SD::CONE:
        ret = {
            {
                .name    = "Half angle",
                .content = M_PI / 4.0,
                .min     = 0.0,
                .max     = M_PI / 2.0,
                .type    = AngleSurfaceParameter,
            },
        };
        break;
    case SD::CYLINDER:
        ret = {
            {
                .name    = "Vertex curvature",
                .content = 1.0,
                .min     = 0.0,
            },
        };
        break;
    case SD::FLAT: break;
    case SD::PARABOLA:
        ret = {
            { .name = "Vertex curvature X", .content = 1.0 },
            { .name = "Vertex curvature Y", .content = 1.0 },
        };
        break;
    case SD::SPHERE:
        ret = {
            { .name = "Vertex curvature", .content = 1.0 },
        };
        break;
    case SD::HYPER:
    case SD::GENERAL_SPENCER_MURTY:
    case SD::TORUS:
    case SD::SURFACE_UNKNOWN: break;
    }

    return ret;
}

std::vector<double>
factory_arguments_for(SD::SurfaceType                  type,
                      QVector<SurfaceParameter> const& records) {
    std::vector<double> ret;

    switch (type) {
    case SD::CONE:
    case SD::CYLINDER:
    case SD::PARABOLA:
    case SD::SPHERE:
        ret.reserve(records.size());
        for (auto const& record : records) {
            ret.push_back(record.content);
        }
        break;
    case SD::FLAT:
    case SD::HYPER:
    case SD::GENERAL_SPENCER_MURTY:
    case SD::TORUS:
    case SD::SURFACE_UNKNOWN: break;
    }

    return ret;
}

double curvature_from_focal_length(double focal_length) {
    if (focal_length == 0.0) return 0.0;
    return 1.0 / (2.0 * focal_length);
}

double focal_length_from_curvature(double curvature) {
    if (curvature == 0.0) return 0.0;
    return 1.0 / (2.0 * curvature);
}

} // namespace

SurfaceParameterModel::SurfaceParameterModel(QObject* parent)
    : StructTableModel(parent),
      m_surface_type_model(new QStringListModel(this)) {
    build_options<SD::SurfaceType>(*m_surface_type_model);

    connect(this,
            &SurfaceParameterModel::dataChanged,
            this,
            &SurfaceParameterModel::updated);

    connect(this,
            &SurfaceParameterModel::dataChanged,
            this,
            [this](auto const&, auto const&, auto const&) {
                if (!database()) return;

                database()->geometry_parameters.try_patch(
                    m_current_group, [this](GeometryComponent& params) {
                        if (params.surface) { write_back(*params.surface); }
                    });
            });

    connect(this,
            &SurfaceParameterModel::surface_kind_changed,
            this,
            &SurfaceParameterModel::updated);

    connect(this,
            &SurfaceParameterModel::surface_kind_changed,
            this,
            &SurfaceParameterModel::surf_changed);
}

void SurfaceParameterModel::set(Database* database, entt::entity group) {
    observe(database);
    m_current_group = group;

    parameters_changed(group);
}

void SurfaceParameterModel::set_new_database_connections(Database* ptr) {
    add_connection(connect(ptr->geometry_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &SurfaceParameterModel::parameters_changed));
}

void SurfaceParameterModel::parameters_changed(entt::entity e) {
    if (!this->database()) return;
    if (m_current_group != e) return;

    auto* params = database()->geometry_parameters.get(m_current_group);
    if (!params) return;

    auto new_surf = SD::SurfaceType::FLAT;

    if (params->surface) {
        new_surf = params->surface->my_type;
        set_from(*params->surface);
    } else {
        set_for(new_surf);
    }

    m_syncing_from_database = true;
    set_surface_kind(QString(magic_enum::enum_name(new_surf).data()));
    m_syncing_from_database = false;
}

void SurfaceParameterModel::set_for(SD::SurfaceType type) {
    this->reset(make_parameters_for(type));
}

void SurfaceParameterModel::set_from(SD::Surface const& surface) {
    auto ret = make_parameters_for(surface.my_type);

    switch (surface.my_type) {
    case SD::CONE: {
        auto const* ptr = dynamic_cast<SD::Cone const*>(&surface);
        if (!ptr) break;
        ret[0].content = ptr->half_angle;
        break;
    }
    case SD::CYLINDER: {
        auto const* ptr = dynamic_cast<SD::Cylinder const*>(&surface);
        if (!ptr) break;
        ret[0].content = ptr->radius == 0.0 ? 0.0 : 1.0 / ptr->radius;
        break;
    }
    case SD::FLAT: break;
    case SD::PARABOLA: {
        auto const* ptr = dynamic_cast<SD::Parabola const*>(&surface);
        if (!ptr) break;
        ret[0].content = curvature_from_focal_length(ptr->focal_length_x);
        ret[1].content = curvature_from_focal_length(ptr->focal_length_y);
        break;
    }
    case SD::SPHERE: {
        auto const* ptr = dynamic_cast<SD::Sphere const*>(&surface);
        if (!ptr) break;
        ret[0].content = ptr->vertex_curv;
        break;
    }
    case SD::HYPER:
    case SD::GENERAL_SPENCER_MURTY:
    case SD::TORUS:
    case SD::SURFACE_UNKNOWN: break;
    }

    this->reset(ret);
}

void SurfaceParameterModel::write_back(SD::Surface& surface) const {
    switch (surface.my_type) {
    case SD::CONE: {
        auto* ptr = dynamic_cast<SD::Cone*>(&surface);
        if (!ptr || m_records.isEmpty()) break;
        ptr->half_angle = m_records[0].content;
        break;
    }
    case SD::CYLINDER: {
        auto* ptr = dynamic_cast<SD::Cylinder*>(&surface);
        if (!ptr || m_records.isEmpty()) break;
        ptr->radius =
            m_records[0].content == 0.0 ? 0.0 : 1.0 / m_records[0].content;
        break;
    }
    case SD::FLAT: break;
    case SD::PARABOLA: {
        auto* ptr = dynamic_cast<SD::Parabola*>(&surface);
        if (!ptr || m_records.size() < 2) break;
        ptr->focal_length_x = focal_length_from_curvature(m_records[0].content);
        ptr->focal_length_y = focal_length_from_curvature(m_records[1].content);
        break;
    }
    case SD::SPHERE: {
        auto* ptr = dynamic_cast<SD::Sphere*>(&surface);
        if (!ptr || m_records.isEmpty()) break;
        ptr->vertex_curv = m_records[0].content;
        break;
    }
    case SD::HYPER:
    case SD::GENERAL_SPENCER_MURTY:
    case SD::TORUS:
    case SD::SURFACE_UNKNOWN: break;
    }
}

void SurfaceParameterModel::make_new_surface(SD::SurfaceType type) {
    if (!database()) return;

    set_for(type);

    auto replacement = SD::make_surface_from_type(
        type, factory_arguments_for(type, m_records));

    if (!replacement) {
        parameters_changed(m_current_group);
        return;
    }

    database()->geometry_parameters.try_patch(
        m_current_group,
        [&](GeometryComponent& params) { params.surface = replacement; });
}

void SurfaceParameterModel::surf_changed() {
    if (m_syncing_from_database) return;

    auto str = surface_kind().toStdString();

    auto new_surf = magic_enum::enum_cast<SD::SurfaceType>(str).value_or(
        SD::SurfaceType::SURFACE_UNKNOWN);

    set_surface_kind(QString(magic_enum::enum_name(new_surf).data()));
    make_new_surface(new_surf);
}

} // namespace db
