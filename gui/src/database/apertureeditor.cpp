#include "apertureeditor.h"

#include <magic_enum/magic_enum.hpp>

#include <QtMath>

namespace db {

namespace {

QVector<ApertureParameter> make_parameters_for(SD::ApertureType type) {
    QVector<ApertureParameter> ret;

    switch (type) {
    case SolTrace::Data::ANNULUS:
        ret = {
            {
                .name    = "Inner radius",
                .content = 0.0,
                .min     = 0.0,
            },
            {
                .name    = "Outer radius",
                .content = 1.0,
                .min     = 0.0,
            },
            {
                .name    = "Arc angle",
                .content = 2 * M_PI,
                .min     = 0.0,
                .max     = 2 * M_PI,
                .type    = AngleApertureParameter,
            },
        };
        break;
    case SolTrace::Data::CIRCLE:
        ret = {
            {
                .name    = "Diameter",
                .content = 1.0,
                .min     = 0.0,
            },
        };
        break;
    case SolTrace::Data::HEXAGON:
        ret = {
            {
                .name    = "Circumscribe diameter",
                .content = 1.0,
                .min     = 0.0,
            },
        };
        break;
    case SolTrace::Data::RECTANGLE:
        ret = {
            {
                .name    = "Lower-left corner (x)",
                .content = -0.5,
            },
            {
                .name    = "Lower-left corner (y)",
                .content = -0.5,
            },
            {
                .name    = "Size (x)",
                .content = 1.0,
                .min     = 0.0,
            },
            {
                .name    = "Size (y)",
                .content = 1.0,
                .min     = 0.0,
            },
        };
        break;
    case SolTrace::Data::EQUILATERAL_TRIANGLE:
        ret = {
            {
                .name    = "Circumscribe diameter",
                .content = 1.0,
                .min     = 0.0,
            },
        };
        break;
    case SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION:
        // ret = {
        //     {
        //         .name    = "X min",
        //         .content = -0.5,
        //     },
        //     {
        //         .name    = "X max",
        //         .content = 0.5,
        //     },
        //     {
        //         .name    = "Y length (x)",
        //         .content = 1.0,
        //         .min     = 0.0,
        //     },
        //     {
        //         .name    = "Y length (y)",
        //         .content = 0.0,
        //         .min     = 0.0,
        //     },
        // };
        break;
    case SolTrace::Data::IRREGULAR_TRIANGLE:
        ret = {
            {
                .name    = "Point 1 (x)",
                .content = 0.0,
            },
            {
                .name    = "Point 1 (y)",
                .content = 1.0,
            },
            {
                .name    = "Point 2 (x)",
                .content = 0.0,
            },
            {
                .name    = "Point 2 (y)",
                .content = 0.0,
            },
            {
                .name    = "Point 3 (x)",
                .content = 1.0,
            },
            {
                .name    = "Point 3 (y)",
                .content = 0.0,
            },
        };
        break;
    case SolTrace::Data::IRREGULAR_QUADRILATERAL:
        ret = {
            {
                .name    = "Point 1 (x)",
                .content = -1.0,
            },
            {
                .name    = "Point 1 (y)",
                .content = -1.0,
            },
            {
                .name    = "Point 2 (x)",
                .content = -1.0,
            },
            {
                .name    = "Point 2 (y)",
                .content = 1.0,
            },
            {
                .name    = "Point 3 (x)",
                .content = 1.0,
            },
            {
                .name    = "Point 3 (y)",
                .content = 1.0,
            },
            {
                .name    = "Point 4 (x)",
                .content = 1.0,
            },
            {
                .name    = "Point 4 (y)",
                .content = -1.0,
            },
        };
        break;
    case SolTrace::Data::APERTURE_UNKNOWN: break;
    }

    return ret;
}

std::vector<double> factory_arguments_for(
    SD::ApertureType type, QVector<ApertureParameter> const& records) {
    if (type == SD::ApertureType::RECTANGLE && records.size() >= 4) {
        return { records[2].content, records[3].content };
    }

    std::vector<double> ret;

    for (auto const& record : records) {
        ret.push_back(record.content);
    }

    return ret;
}

} // namespace

std::span<SD::ApertureType const>
ApertureParameterModel::valid_apertures_for_surf(SD::SurfaceType st) {

    static constexpr std::array<SD::ApertureType, 8> all_valid = {
        SolTrace::Data::ANNULUS,
        SolTrace::Data::CIRCLE,
        SolTrace::Data::HEXAGON,
        SolTrace::Data::RECTANGLE,
        SolTrace::Data::EQUILATERAL_TRIANGLE,
        SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION,
        SolTrace::Data::IRREGULAR_TRIANGLE,
        SolTrace::Data::IRREGULAR_QUADRILATERAL,
    };

    switch (st) {
    case SolTrace::Data::CYLINDER: {
        static constexpr std::array<SD::ApertureType, 1> valid = {
            SolTrace::Data::RECTANGLE,
        };
        return valid;
    }
    case SolTrace::Data::HYPER: {
        static constexpr std::array<SD::ApertureType, 5> valid = {
            SolTrace::Data::ANNULUS,
            SolTrace::Data::CIRCLE,
            SolTrace::Data::HEXAGON,
            SolTrace::Data::RECTANGLE,
            SolTrace::Data::EQUILATERAL_TRIANGLE,
        };
        return valid;
    }
    case SolTrace::Data::CONE:
    case SolTrace::Data::FLAT:
    case SolTrace::Data::PARABOLA:
    case SolTrace::Data::SPHERE:
    case SolTrace::Data::GENERAL_SPENCER_MURTY:
    case SolTrace::Data::TORUS: return all_valid;
    case SolTrace::Data::SURFACE_UNKNOWN: return {};
    }

    return {};
}

ApertureParameterModel::ApertureParameterModel(QObject* parent)
    : StructTableModel(parent),
      m_aperture_type_model(new QStringListModel(this)) {

    build_options<SD::ApertureType>(*m_aperture_type_model);

    connect(this,
            &ApertureParameterModel::dataChanged,
            this,
            &ApertureParameterModel::updated);

    connect(this,
            &ApertureParameterModel::dataChanged,
            this,
            [this](auto const&, auto const&, auto const&) {
                if (!database()) return;

                database()->geometry_parameters.try_patch(
                    m_current_group, [this](GeometryComponent& params) {
                        if (params.aperture) {
                            write_back(*params.aperture);
                        }
                    });
            });

    connect(this,
            &ApertureParameterModel::aperture_kind_changed,
            this,
            &ApertureParameterModel::updated);

    connect(this,
            &ApertureParameterModel::aperture_kind_changed,
            this,
            &ApertureParameterModel::apt_changed);
}

void ApertureParameterModel::set(Database* database, entt::entity group) {
    observe(database);
    m_current_group = group;

    parameters_changed(group);
}

void ApertureParameterModel::set_new_database_connections(Database* ptr) {
    add_connection(connect(ptr->geometry_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &ApertureParameterModel::parameters_changed));
}

void ApertureParameterModel::parameters_changed(entt::entity e) {
    if (!this->database()) return;
    if (m_current_group != e) return;

    auto* params = database()->geometry_parameters.get(m_current_group);
    if (!params) return;

    auto new_apt = SD::ApertureType::RECTANGLE;

    if (params->aperture) {
        new_apt = params->aperture->my_type;
        set_from(*params->aperture);
    } else {
        set_for(new_apt);
    }

    m_syncing_from_database = true;
    set_aperture_kind(QString(magic_enum::enum_name(new_apt).data()));
    m_syncing_from_database = false;
}

void ApertureParameterModel::set_for(SD::ApertureType type) {
    this->reset(make_parameters_for(type));
}

void ApertureParameterModel::set_from(SD::Aperture const& aperture) {
    auto ret = make_parameters_for(aperture.my_type);

    switch (aperture.my_type) {
    case SolTrace::Data::ANNULUS: {
        auto const* ptr = dynamic_cast<SD::Annulus const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->inner_radius;
        ret[1].content = ptr->outer_radius;
        ret[2].content = ptr->arc_angle;
        break;
    }
    case SolTrace::Data::CIRCLE: {
        auto const* ptr = dynamic_cast<SD::Circle const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->diameter;
        break;
    }
    case SolTrace::Data::HEXAGON: {
        auto const* ptr = dynamic_cast<SD::Hexagon const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->circumscribe_diameter;
        break;
    }
    case SolTrace::Data::RECTANGLE: {
        auto const* ptr = dynamic_cast<SD::Rectangle const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->x_coord();
        ret[1].content = ptr->y_coord();

        ret[2].content = ptr->x_length();
        ret[3].content = ptr->y_length();
        break;
    }
    case SolTrace::Data::EQUILATERAL_TRIANGLE: {
        auto const* ptr =
            dynamic_cast<SD::EquilateralTriangle const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->circumscribe_diameter;
        break;
    }
    case SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION:
        // TODO: Wait for class to be implemented
        break;
    case SolTrace::Data::IRREGULAR_TRIANGLE: {
        auto const* ptr = dynamic_cast<SD::IrregularTriangle const*>(&aperture);
        if (!ptr) break;

        ret[0].content = ptr->x1;
        ret[1].content = ptr->y1;

        ret[2].content = ptr->x2;
        ret[3].content = ptr->y2;

        ret[4].content = ptr->x3;
        ret[5].content = ptr->y3;

        break;
    }
    case SolTrace::Data::IRREGULAR_QUADRILATERAL: {
        auto const* ptr =
            dynamic_cast<SD::IrregularQuadrilateral const*>(&aperture);
        if (!ptr) break;
        ret[0].content = ptr->x1;
        ret[1].content = ptr->y1;

        ret[2].content = ptr->x2;
        ret[3].content = ptr->y2;

        ret[4].content = ptr->x3;
        ret[5].content = ptr->y3;

        ret[6].content = ptr->x4;
        ret[7].content = ptr->y4;
        break;
    }
    case SolTrace::Data::APERTURE_UNKNOWN: break;
    }

    this->reset(ret);
}

void ApertureParameterModel::write_back(SD::Aperture& aperture) const {
    switch (aperture.my_type) {
    case SolTrace::Data::ANNULUS: {
        auto* ptr = dynamic_cast<SD::Annulus*>(&aperture);
        if (!ptr || m_records.size() < 3) break;
        ptr->inner_radius = m_records[0].content;
        ptr->outer_radius = m_records[1].content;
        ptr->arc_angle    = m_records[2].content;
        break;
    }
    case SolTrace::Data::CIRCLE: {
        auto* ptr = dynamic_cast<SD::Circle*>(&aperture);
        if (!ptr || m_records.isEmpty()) break;
        ptr->diameter = m_records[0].content;
        break;
    }
    case SolTrace::Data::HEXAGON: {
        auto* ptr = dynamic_cast<SD::Hexagon*>(&aperture);
        if (!ptr || m_records.isEmpty()) break;
        ptr->circumscribe_diameter = m_records[0].content;
        break;
    }
    case SolTrace::Data::RECTANGLE: {
        auto* ptr = dynamic_cast<SD::Rectangle*>(&aperture);
        if (!ptr || m_records.size() < 4) break;
        ptr->set_x_coord(m_records[0].content);
        ptr->set_y_coord(m_records[1].content);
        ptr->set_x_length(m_records[2].content);
        ptr->set_y_length(m_records[3].content);
        break;
    }
    case SolTrace::Data::EQUILATERAL_TRIANGLE: {
        auto* ptr = dynamic_cast<SD::EquilateralTriangle*>(&aperture);
        if (!ptr || m_records.isEmpty()) break;
        ptr->circumscribe_diameter = m_records[0].content;
        break;
    }
    case SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION:
        break;
    case SolTrace::Data::IRREGULAR_TRIANGLE: {
        auto* ptr = dynamic_cast<SD::IrregularTriangle*>(&aperture);
        if (!ptr || m_records.size() < 6) break;
        ptr->x1 = m_records[0].content;
        ptr->y1 = m_records[1].content;
        ptr->x2 = m_records[2].content;
        ptr->y2 = m_records[3].content;
        ptr->x3 = m_records[4].content;
        ptr->y3 = m_records[5].content;
        break;
    }
    case SolTrace::Data::IRREGULAR_QUADRILATERAL: {
        auto* ptr = dynamic_cast<SD::IrregularQuadrilateral*>(&aperture);
        if (!ptr || m_records.size() < 8) break;
        ptr->x1 = m_records[0].content;
        ptr->y1 = m_records[1].content;
        ptr->x2 = m_records[2].content;
        ptr->y2 = m_records[3].content;
        ptr->x3 = m_records[4].content;
        ptr->y3 = m_records[5].content;
        ptr->x4 = m_records[6].content;
        ptr->y4 = m_records[7].content;
        break;
    }
    case SolTrace::Data::APERTURE_UNKNOWN:
        break;
    }
}

void ApertureParameterModel::make_new_aperture(SD::ApertureType type) {
    if (!database()) return;

    set_for(type);

    auto replacement = SD::Aperture::make_aperture_from_type(
        type, factory_arguments_for(type, m_records));

    if (!replacement) {
        parameters_changed(m_current_group);
        return;
    }

    database()->geometry_parameters.try_patch(
        m_current_group,
        [&](GeometryComponent& params) { params.aperture = replacement; });
}

void ApertureParameterModel::apt_changed() {
    if (m_syncing_from_database) return;

    auto str = aperture_kind().toStdString();

    auto new_apt = magic_enum::enum_cast<SD::ApertureType>(str).value_or(
        SD::ApertureType::APERTURE_UNKNOWN);

    set_aperture_kind(QString(magic_enum::enum_name(new_apt).data()));
    make_new_aperture(new_apt);
}

} // namespace db
