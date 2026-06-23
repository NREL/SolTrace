#include "geometryeditor.h"

#include "database/apertureeditor.h"
#include "database/components.h"
#include "database/surface.h"

#include <magic_enum/magic_enum.hpp>

#include <algorithm>
#include <cmath>

namespace db {

SurfaceGeometry::SurfaceGeometry() {
    connect(this,
            &SurfaceGeometry::quality_changed,
            this,
            &SurfaceGeometry::rebuild_geometry);
    connect(this,
            &SurfaceGeometry::add_thickness_changed,
            this,
            &SurfaceGeometry::rebuild_geometry);
    connect(this,
            &SurfaceGeometry::thickness_changed,
            this,
            &SurfaceGeometry::rebuild_geometry);

    rebuild_geometry();
}

void SurfaceGeometry::set_new_database_connections(Database const* ptr) {

    // careful here
    auto* mptr = const_cast<Database*>(ptr);

    add_connection(connect(mptr->geometry_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &SurfaceGeometry::parameters_changed));
}

void SurfaceGeometry::set(Database const* ptr, entt::entity group) {
    observe(ptr);
    m_current_group = group;

    rebuild_geometry();
}

void SurfaceGeometry::parameters_changed(entt::entity group) {
    if (group == m_current_group) rebuild_geometry();
}

SurfaceGenerationOptions SurfaceGeometry::surface_generation_options() const {
    SurfaceGenerationOptions options;

    auto scale = [this](uint32_t value) {
        switch (quality()) {
        case Quality::Low:
            return std::max<uint32_t>(1, value / 2);
        case Quality::Normal:
            return value;
        case Quality::High:
            return value * 2;
        }

        return value;
    };

    options.height_field_resolution = {
        scale(options.height_field_resolution.x),
        scale(options.height_field_resolution.y),
    };
    options.radial_subdivisions           = scale(options.radial_subdivisions);
    options.perimeter_subdivisions        = scale(options.perimeter_subdivisions);
    options.cylinder_angular_subdivisions = scale(options.cylinder_angular_subdivisions);
    options.cylinder_length_subdivisions  = scale(options.cylinder_length_subdivisions);
    options.add_thickness                 = add_thickness();
    options.thickness                     = thickness();

    return options;
}

void SurfaceGeometry::rebuild_geometry() {
    clear();

    if (!database()) {
        qDebug() << Q_FUNC_INFO << "no db";
        return;
    }

    auto ptr = database()->geometry_parameters.get(m_current_group);

    if (!ptr) {
        qDebug() << Q_FUNC_INFO << "no geometry";
        return;
    }

    if (!ptr->surface || !ptr->aperture) {
        qDebug() << Q_FUNC_INFO << "no surf or apt";
        update();
        return;
    }

    auto surface  = ptr->surface;
    auto aperture = ptr->aperture;

    auto mesh = generate_surface(surface, aperture, surface_generation_options());
    if (!mesh || mesh->vertex.empty() || mesh->triangles.empty()) {
        qWarning() << Q_FUNC_INFO
                   << "Geometry is empty, or unable to be generated";
        set_vertex_count(0);
        update();
        return;
    }

    constexpr float max_float = std::numeric_limits<float>::max();

    glm::vec3 bounds_min(max_float);
    glm::vec3 bounds_max(-max_float);

    for (auto const& p : std::as_const(mesh->vertex)) {
        bounds_min = glm::min(bounds_min, p.position);
        bounds_max = glm::max(bounds_max, p.position);
    }

    auto extent     = bounds_max - bounds_min;
    auto max_extent = std::max({ extent.x, extent.y, extent.z, 1.0f });
    auto padding    = max_extent * 1.0e-4f;

    for (int axis = 0; axis < 3; ++axis) {
        if (bounds_min[axis] == bounds_max[axis]) {
            bounds_min[axis] -= padding;
            bounds_max[axis] += padding;
        }
    }

    auto indexBuffer =
        QByteArray(reinterpret_cast<const char*>(mesh->triangles.data()),
                   mesh->triangles.size() * sizeof(glm::uvec3));
    auto vertexBuffer =
        QByteArray(reinterpret_cast<const char*>(mesh->vertex.data()),
                   mesh->vertex.size() * sizeof(Vertex));

    addAttribute(QQuick3DGeometry::Attribute::PositionSemantic,
                 offsetof(Vertex, position),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::NormalSemantic,
                 offsetof(Vertex, normal),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::TexCoord0Semantic,
                 offsetof(Vertex, uv),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::IndexSemantic,
                 0,
                 QQuick3DGeometry::Attribute::ComponentType::U32Type);

    auto bb = BoundingBox {
        .min = QVector3D(bounds_min.x, bounds_min.y, bounds_min.z),
        .max = QVector3D(bounds_max.x, bounds_max.y, bounds_max.z),
    };

    setStride(sizeof(Vertex));
    setVertexData(vertexBuffer);
    setIndexData(indexBuffer);
    setBounds(bb.min, bb.max);
    setPrimitiveType(QQuick3DGeometry::PrimitiveType::Triangles);

    set_bounding_box(bb);

    set_vertex_count(mesh->vertex.count());

    // qDebug() << Q_FUNC_INFO << entt::to_integral(m_current_group)
    //          << mesh->triangles.size() << mesh->vertex.size();
    // qDebug() << Q_FUNC_INFO << bb.min << bb.max;
    //   qDebug() << verts;

    update();
}

void SurfaceGeometry::debug() {
    qDebug() << this->bounding_box().min << this->bounding_box().max;
}

GeometryEditor::GeometryEditor(QObject* parent)
    : QObject { parent },
      m_surface_geometry(new SurfaceGeometry()),
      m_aperture_parameter_model(new ApertureParameterModel(this)),
      m_surface_parameter_model(new SurfaceParameterModel(this)),
      m_geometry_error_model(new QStringListModel(this))

{
    m_surface_geometry->setParent(this);
    m_surface_geometry->set_quality(SurfaceGeometry::Quality::High);
    m_surface_geometry->set_add_thickness(true);
    m_surface_geometry->set_thickness(.01);

    connect(m_surface_parameter_model,
            &SurfaceParameterModel::updated,
            this,
            &GeometryEditor::updated);

    connect(m_aperture_parameter_model,
            &ApertureParameterModel::updated,
            this,
            &GeometryEditor::updated);
}

GeometryEditor::~GeometryEditor() {
}

void GeometryEditor::set_new_database_connections(Database* ptr) {
    add_connection(connect(ptr->geometry_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &GeometryEditor::geometry_parameters_changed));
}

void GeometryEditor::geometry_parameters_changed(entt::entity group) {
    if (group == m_current_group) recompute_geometry_errors();
}

void GeometryEditor::recompute_geometry_errors() {
    QStringList errors;

    auto* params = database() ? database()->geometry_parameters.get(m_current_group)
                              : nullptr;
    if (!params) {
        m_geometry_error_model->setStringList(errors);
        return;
    }

    if (!params->surface || !params->aperture) {
        errors.push_back("Geometry is missing a surface or aperture.");
        m_geometry_error_model->setStringList(errors);
        return;
    }

    auto valid_apertures =
        ApertureParameterModel::valid_apertures_for_surf(params->surface->my_type);
    auto is_valid = std::find(valid_apertures.begin(),
                              valid_apertures.end(),
                              params->aperture->my_type) != valid_apertures.end();

    if (!is_valid) {
        auto surface_name =
            QString(magic_enum::enum_name(params->surface->my_type).data());
        auto aperture_name =
            QString(magic_enum::enum_name(params->aperture->my_type).data());

        errors.push_back(
            QString("%1 is not a valid aperture for %2 surfaces.")
                .arg(aperture_name, surface_name));
    }

    if (params->surface->my_type == SD::CYLINDER &&
        params->aperture->my_type == SD::RECTANGLE) {

        auto const* cylinder =
            dynamic_cast<SD::Cylinder const*>(params->surface.get());
        auto const* rect =
            dynamic_cast<SD::Rectangle const*>(params->aperture.get());

        if (cylinder && rect) {
            auto close = [](double a, double b) {
                return std::abs(a - b) <= 1.0e-8;
            };

            const double diameter = 2.0 * cylinder->radius;

            if (cylinder->radius <= 0.0) {
                errors.push_back("Cylinder radius must be positive.");
            }

            if (!close(rect->x_length(), diameter)) {
                errors.push_back(
                    "Cylinder rectangle width must equal the cylinder diameter.");
            }

            if (!close(rect->x_coord(), -cylinder->radius)) {
                errors.push_back(
                    "Cylinder rectangle X origin must be the negative radius.");
            }
        }
    }

    m_geometry_error_model->setStringList(errors);
}

void GeometryEditor::set(Database* database, entt::entity group) {
    observe(database);
    m_current_group = group;
    m_surface_geometry->set(database, group);
    m_aperture_parameter_model->set(database, group);
    m_surface_parameter_model->set(database, group);
    recompute_geometry_errors();
}

} // namespace db
