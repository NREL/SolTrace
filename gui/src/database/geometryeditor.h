#pragma once

#include <QObject>
#include <QStringListModel>
#include <QtGui/qvectornd.h>
#include <QtQuick3D/qquick3dgeometry.h>

#include "database/apertureeditor.h"
#include "database/components.h"
#include "database/database.h"
#include "database/database_observer.h"
#include "database/models/surface_parameter_model.h"
#include "utilities/qt_helpers.h"

#include "aperture.hpp"
#include "surface.hpp"

namespace SD = SolTrace::Data;


namespace db {

struct MaterialComponent;
struct SurfaceGenerationOptions;

/// Axis-aligned bounds for generated surface geometry.
class BoundingBox {
    Q_GADGET
    Q_PROPERTY(QVector3D min MEMBER min)
    Q_PROPERTY(QVector3D max MEMBER max)

public:
    QVector3D min;
    QVector3D max;

    bool operator==(BoundingBox const&) const = default;
};


/// Surface geometry visualization for one geometry group.
///
/// Rebuilds Quick3D geometry buffers when surface or aperture parameters
/// change in the observed database.
class SurfaceGeometry : public QQuick3DGeometry, public ConstDatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group = entt::null;

    void set_new_database_connections(Database const* ptr) override;
    SurfaceGenerationOptions surface_generation_options() const;

private slots:
    void parameters_changed(entt::entity);
    void rebuild_geometry();

    Q_READONLY_PROPERTY(unsigned, vertex_count);

public:
    enum class Quality { Low, Normal, High };
    Q_ENUM(Quality)

    SurfaceGeometry();

    /// The quality of the surface geometry. Rough control on surface
    /// subdivision
    Q_WRITABLE_PROPERTY(Quality, quality, Quality::Normal)

    /// Add thickness to the surface (ie, doubleside with edges)
    Q_WRITABLE_PROPERTY(bool, add_thickness, false)

    /// If thickness is added, how much in world units
    Q_WRITABLE_PROPERTY(double, thickness, 0.01)

    /// Subdiv control. TODO: Clarify with above
    Q_WRITABLE_PROPERTY(unsigned, subdivision_scale, 2)

    /// Bounding box of the surface geometry item
    Q_READONLY_PROPERTY(BoundingBox, bounding_box)

    /// Observe database geometry group and rebuild the generated geometry.
    void set(Database const*, entt::entity group);

public:
    /// Print debugging information about the generated geometry.
    void debug();
};


/// QML-facing editor for a geometry group.
class GeometryEditor : public QObject, public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group = entt::null;

    void set_new_database_connections(Database* ptr) override;

    QOBJECT_WRITABLE_PROPERTY(SurfaceGeometry, surface_geometry);

private:
    QOBJECT_READONLY_PROPERTY(ApertureParameterModel, aperture_parameter_model);
    QOBJECT_READONLY_PROPERTY(SurfaceParameterModel, surface_parameter_model);
    QOBJECT_READONLY_PROPERTY(QStringListModel, geometry_error_model);

private slots:
    void geometry_parameters_changed(entt::entity);
    void recompute_geometry_errors();

public:
    explicit GeometryEditor(QObject* parent = nullptr);
    ~GeometryEditor() override;

    /// Observe database geometry group and synchronize parameter models.
    void set(Database*, entt::entity group);

signals:
    void updated();
};

} // namespace db

Q_DECLARE_METATYPE(db::BoundingBox)
