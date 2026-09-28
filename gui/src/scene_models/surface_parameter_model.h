#pragma once

#include <QObject>
#include <QStringListModel>
#include <limits>
#include <vector>

#include "surface.hpp"

#include "data/database.h"
#include "data/database_observer.h"
#include "support/qt_helpers.h"
#include "support/structmodel.h"

namespace SD = SolTrace::Data;

namespace SolTrace::GUI::Data {

/// UI presentation hint for a surface parameter.
enum SurfaceParameterType {
    GenericSurfaceParameter = 0,
    AngleSurfaceParameter   = 1,
};

/// One editable surface parameter exposed to QML.
struct SurfaceParameter {
    QString name;
    double  content = 0.0;
    double  min     = std::numeric_limits<double>::lowest();
    double  max     = std::numeric_limits<double>::max();
    int     type    = GenericSurfaceParameter;


    RECORD_META(SurfaceParameter,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RW(content),
                SM_EXPOSE_RO(min),
                SM_EXPOSE_RO(max),
                SM_EXPOSE_RO(type));
};


/// Model that adapts SolTrace surface parameters to editable QML rows.
class SurfaceParameterModel
    : public Support::StructTableModel<SurfaceParameter>,
      public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group         = entt::null;
    bool         m_syncing_from_database = false;

    void set_new_database_connections(Database* ptr) override;

    Q_WRITABLE_PROPERTY(QString, surface_kind, "SURFACE_UNKNOWN");
    QOBJECT_READONLY_PROPERTY(QStringListModel, surface_type_model);

    void make_new_surface(SD::SurfaceType);

public:
    explicit SurfaceParameterModel(QObject* parent = nullptr);

    /// Observe a database geometry group.
    void set(Database*, entt::entity group);

    /// Replace the current surface with a new surface of type.
    void set_for(SD::SurfaceType);

    /// Populate rows from an existing SolTrace surface.
    void set_from(SD::Surface const&);

    /// Write edited rows back into a SolTrace surface object.
    void write_back(SD::Surface&) const;

private slots:
    void parameters_changed(entt::entity);
    void surf_changed();

signals:
    void updated();
};

} // namespace SolTrace::GUI::Data
