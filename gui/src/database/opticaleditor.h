#pragma once

#include "database/database.h"
#include "simulation_data_export.hpp"

#include <QObject>
#include <QtCore/qstringlistmodel.h>
#include <entt/entity/fwd.hpp>

namespace db {

class OpticalPropertiesObject : public QObject, public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group = entt::null;
    bool         m_back          = false;

    SD::OpticalSide m_side = SD::OpticalSide::Front;

    SolTrace::Data::OpticalPropertySet const* get_properties() const;

    template <class F>
    void patch_properties(F&& f) {
        if (!database()) return;

        database()->material_parameters.try_patch(
            m_current_group,
            [this, &f](MaterialComponent& c) { f(c.optics, m_side); });
    }

private:
    // TODO: interaction is now unified front and back.

    Q_PROPERTY(QString interaction_type READ interaction_type WRITE
                   set_interaction_type NOTIFY interaction_type_changed)

    Q_PROPERTY(
        QString error_distribution_type READ error_distribution_type WRITE
            set_error_distribution_type NOTIFY error_distribution_type_changed)

    Q_PROPERTY(double transmissivity READ transmissivity WRITE
                   set_transmissivity NOTIFY transmissivity_changed)

    Q_PROPERTY(double reflectivity READ reflectivity WRITE set_reflectivity
                   NOTIFY reflectivity_changed)

    Q_PROPERTY(double slope_error READ slope_error WRITE set_slope_error NOTIFY
                   slope_error_changed)

    Q_PROPERTY(double specularity_error READ specularity_error WRITE
                   set_specularity_error NOTIFY specularity_error_changed)

    Q_PROPERTY(
        double refraction_index_front READ refraction_index_front WRITE
            set_refraction_index_front NOTIFY refraction_index_front_changed)

    Q_PROPERTY(
        double refraction_index_back READ refraction_index_back WRITE
            set_refraction_index_back NOTIFY refraction_index_back_changed)

    // UX helpers
    Q_PROPERTY(double transmissivity_min READ transmissivity_min CONSTANT)
    Q_PROPERTY(double transmissivity_max READ transmissivity_max CONSTANT)
    Q_PROPERTY(double reflectivity_min READ reflectivity_min CONSTANT)
    Q_PROPERTY(double reflectivity_max READ reflectivity_max CONSTANT)

    void set_new_database_connections(Database* ptr) override;

    void trigger_all_changed();

private slots:
    void parameters_changed(entt::entity);

public:
    explicit OpticalPropertiesObject(bool back, QObject* parent = nullptr);
    ~OpticalPropertiesObject() override;

    void set(Database*, entt::entity group);

public:
    QString interaction_type() const;
    QString error_distribution_type() const;

    double transmissivity() const;
    double reflectivity() const;
    double slope_error() const;
    double specularity_error() const;
    double refraction_index_front() const;
    double refraction_index_back() const;

public slots:
    void set_interaction_type(QString v);
    void set_error_distribution_type(QString v);

    void set_transmissivity(double v);
    void set_reflectivity(double v);
    void set_slope_error(double v);
    void set_specularity_error(double v);
    void set_refraction_index_front(double v);
    void set_refraction_index_back(double v);

    // presets
    void set_ideal_absorption();
    void set_ideal_reflection();
    void set_ideal_transmission();
    void set_ideal_transmission_with_indices(double n_front, double n_back);

signals:
    void interaction_type_changed();
    void error_distribution_type_changed();

    void transmissivity_changed();
    void reflectivity_changed();
    void slope_error_changed();
    void specularity_error_changed();
    void refraction_index_front_changed();
    void refraction_index_back_changed();

private:
    double transmissivity_min() const { return 0.0; }
    double transmissivity_max() const { return 1.0; }
    double reflectivity_min() const { return 0.0; }
    double reflectivity_max() const { return 1.0; }
};

} // namespace db
