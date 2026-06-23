#include "opticaleditor.h"

#include <magic_enum/magic_enum.hpp>

namespace db {

OpticalPropertiesObject::OpticalPropertiesObject(bool back, QObject* parent)
    : QObject(parent),
      m_current_group(entt::null),
      m_back(back),
      m_side(m_back ? SD::OpticalSide::Back : SD::OpticalSide::Front) { }

OpticalPropertiesObject::~OpticalPropertiesObject() = default;

void OpticalPropertiesObject::trigger_all_changed() {
    qDebug() << Q_FUNC_INFO << m_back;
    emit interaction_type_changed();
    emit error_distribution_type_changed();
    emit transmissivity_changed();
    emit reflectivity_changed();
    emit slope_error_changed();
    emit specularity_error_changed();
    emit refraction_index_front_changed();
    emit refraction_index_back_changed();
}

void OpticalPropertiesObject::set(Database* db, entt::entity group) {
    observe(db);

    if (!db or !db->material_parameters.get(group)) {
        m_current_group = entt::null;
        trigger_all_changed();
        return;
    }

    m_current_group = group;

    trigger_all_changed();
}

void OpticalPropertiesObject::set_new_database_connections(Database* ptr) {
    add_connection(connect(ptr->material_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &OpticalPropertiesObject::parameters_changed));
}

void OpticalPropertiesObject::parameters_changed(entt::entity e) {
    if (this->m_current_group != e) return;

    trigger_all_changed();
}

SolTrace::Data::OpticalPropertySet const*
OpticalPropertiesObject::get_properties() const {
    if (!database()) return nullptr;

    auto* ptr = database()->material_parameters.get(m_current_group);

    if (!ptr) return nullptr;

    return &ptr->optics;
}

QString OpticalPropertiesObject::interaction_type() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return QString(magic_enum::enum_name(ptr->get_interaction_type()).data());
}

QString OpticalPropertiesObject::error_distribution_type() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return QString(
        magic_enum::enum_name(ptr->get_error_distribution(m_side)).data());
}

double OpticalPropertiesObject::transmissivity() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return ptr->get_transmissivity(m_side);
}

double OpticalPropertiesObject::reflectivity() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return ptr->get_reflectivity(m_side);
}

double OpticalPropertiesObject::slope_error() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return ptr->get_slope_error(m_side);
}

double OpticalPropertiesObject::specularity_error() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    return ptr->get_specularity_error(m_side);
}

double OpticalPropertiesObject::refraction_index_front() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    double front;
    double back;

    ptr->get_refraction_indices(front, back);

    return front;
}

double OpticalPropertiesObject::refraction_index_back() const {
    auto* ptr = get_properties();

    if (!ptr) return {};

    double front;
    double back;

    ptr->get_refraction_indices(front, back);

    return back;
}


void OpticalPropertiesObject::set_interaction_type(QString v) {
    patch_properties([this, &v](SD::OpticalPropertySet& prop, auto side) {
        auto string = v.toStdString();

        auto value =
            magic_enum::enum_cast<SD::InteractionType>(string).value_or(
                SD::InteractionType::REFLECTION);

        prop.set_interaction_type(value);
        emit interaction_type_changed();
    });
}

void OpticalPropertiesObject::set_error_distribution_type(QString v) {
    patch_properties([this, &v](SD::OpticalPropertySet& prop, auto side) {
        auto string = v.toStdString();

        auto value =
            magic_enum::enum_cast<SD::DistributionType>(string).value_or(
                SD::DistributionType::DIFFUSE);

        auto slope_error = prop.get_slope_error(m_side);
        auto spec_error  = prop.get_specularity_error(m_side);

        prop.set_errors(m_side, value, slope_error, spec_error);

        emit error_distribution_type_changed();
    });
}

#define SETTER(MEMBER)                                                         \
    patch_properties([this, &v](SD::OpticalPropertySet& prop, auto side) {     \
        prop.set_##MEMBER(m_side, v);                                          \
        emit MEMBER##_changed();                                               \
    });

void OpticalPropertiesObject::set_transmissivity(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        prop.set_transmissivity(m_side, v);
        transmissivity_changed();
    });
}

void OpticalPropertiesObject::set_reflectivity(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        prop.set_reflectivity(m_side, v);
        reflectivity_changed();
    });
}

void OpticalPropertiesObject::set_slope_error(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        auto type       = prop.get_error_distribution(m_side);
        auto spec_error = prop.get_specularity_error(m_side);

        prop.set_errors(m_side, type, v, spec_error);

        slope_error_changed();
    });
}

void OpticalPropertiesObject::set_specularity_error(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        auto type        = prop.get_error_distribution(m_side);
        auto slope_error = prop.get_slope_error(m_side);

        prop.set_errors(m_side, type, slope_error, v);

        specularity_error_changed();
    });
}

void OpticalPropertiesObject::set_refraction_index_front(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        double front;
        double back;

        prop.get_refraction_indices(front, back);

        prop.set_refraction_indices(v, back);
        emit refraction_index_front_changed();
    });
}

void OpticalPropertiesObject::set_refraction_index_back(double v) {
    patch_properties([this, &v](SD ::OpticalPropertySet& prop, auto side) {
        double front;
        double back;

        prop.get_refraction_indices(front, back);

        prop.set_refraction_indices(front, v);
        emit refraction_index_back_changed();
    });
}

void OpticalPropertiesObject::set_ideal_absorption() {
    patch_properties([this](SD::OpticalPropertySet& prop, auto side) {
        prop.set_ideal_absorption(m_side);
        trigger_all_changed();
    });
}

void OpticalPropertiesObject::set_ideal_reflection() {
    patch_properties([this](SD::OpticalPropertySet& prop, auto side) {
        prop.set_ideal_reflection(m_side);
        trigger_all_changed();
    });
}

void OpticalPropertiesObject::set_ideal_transmission() {
    patch_properties([this](SD::OpticalPropertySet& prop, auto side) {
        prop.set_ideal_transmission();
        trigger_all_changed();
    });
}

void OpticalPropertiesObject::set_ideal_transmission_with_indices(
    double n_front,
    double n_back) {

    patch_properties(
        [n_front, n_back, this](SD::OpticalPropertySet& prop, auto side) {
            prop.set_ideal_transmission(n_front, n_back);
            trigger_all_changed();
        });
}

} // namespace db
