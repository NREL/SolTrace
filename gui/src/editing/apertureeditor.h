#pragma once

#include <QObject>
#include <QPointF>
#include <QStringListModel>
#include <limits>

#include "aperture.hpp"

#include "data/database.h"
#include "data/database_observer.h"
#include "support/qt_helpers.h"
#include "support/structmodel.h"

namespace SD = SolTrace::Data;

namespace SolTrace::GUI::Data {

/// The type of the aperture parameter
enum ApertureParameterType {
    GenericApertureParameter = 0,
    AngleApertureParameter   = 1,
};

/// An aperture parameter
///
/// TODO: make friendier for coordinates, etc.
struct ApertureParameter {
    QString name;
    double  content = 0.0;
    double  min     = std::numeric_limits<double>::lowest();
    double  max     = std::numeric_limits<double>::max();
    int     type    = GenericApertureParameter;

    RECORD_META(ApertureParameter,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RW(content),
                SM_EXPOSE_RO(min),
                SM_EXPOSE_RO(max),
                SM_EXPOSE_RO(type));
};

/// Model to present aperture parameters to QML
class ApertureParameterModel
    : public Support::StructTableModel<ApertureParameter>,
      public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group         = entt::null;
    bool         m_syncing_from_database = false;

    void set_new_database_connections(Database* ptr) override;

    Q_WRITABLE_PROPERTY(QString, aperture_kind, "APERTURE_UNKNOWN");
    QOBJECT_READONLY_PROPERTY(QStringListModel, aperture_type_model);

    void make_new_aperture(SD::ApertureType);

public:
    /// Get valid aperture types for a surface
    static std::span<SolTrace::Data::ApertureType const>
        valid_apertures_for_surf(SD::SurfaceType);

    explicit ApertureParameterModel(QObject* parent);

    /// Set the database/group to edit
    void set(Database*, entt::entity group);

    /// Update the aperture type (i.e. a new aperture of this type).
    void set_for(SD::ApertureType);

    /// Set from an existing SD aperture
    void set_from(SD::Aperture const&);

    /// Update the SD aperture from this model
    void write_back(SD::Aperture&) const;

private slots:
    void parameters_changed(entt::entity);
    void apt_changed();

signals:
    void updated();
};


} // namespace SolTrace::GUI::Data
