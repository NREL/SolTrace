#pragma once

#include <QObject>
#include <QPointF>
#include <QStringListModel>
#include <limits>
#include <vector>

#include "aperture.hpp"

#include "database/database.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

namespace SD = SolTrace::Data;

namespace db {

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


class ApertureParameterModel : public StructTableModel<ApertureParameter>,
                               public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group = entt::null;
    bool         m_syncing_from_database = false;

    void set_new_database_connections(Database* ptr) override;

    Q_WRITABLE_PROPERTY(QString, aperture_kind, "APERTURE_UNKNOWN");
    QOBJECT_READONLY_PROPERTY(QStringListModel, aperture_type_model);

    void make_new_aperture(SD::ApertureType);

public:
    static std::span<SolTrace::Data::ApertureType const>
        valid_apertures_for_surf(SD::SurfaceType);
    explicit ApertureParameterModel(QObject* parent);

    void set(Database*, entt::entity group);

    void set_for(SD::ApertureType);
    void set_from(SD::Aperture const&);
    void write_back(SD::Aperture&) const;

private slots:
    void parameters_changed(entt::entity);
    void apt_changed();

signals:
    void updated();
};


} // namespace db
