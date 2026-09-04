#include "modules/sun/sun_shape.h"

#include <QFile>
#include <QPointer>
#include <QStringView>
#include <QTextStream>

#include <algorithm>
#include <cmath>

namespace SolTrace::GUI::App {

SunShape::SunShape(QObject* parent)
    : QObject(parent),
      m_generated_distribution(new SunShapeModel(this)),
      m_custom_distribution(new SunShapeModel(this)) {

    // SunShape::shape_changed() -> SunShape::update_current_distribution()
    connect(this,
            &SunShape::shape_changed,
            this,
            &SunShape::update_current_distribution);

    // SunShape::*_changed() -> SunShape::regenerate()
    connect(this, &SunShape::shape_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::sigma_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::half_width_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::csr_changed, this, &SunShape::regenerate);
    connect(m_custom_distribution,
            &SunShapeModel::changed,
            this,
            &SunShape::update_x_axis);

    // SunShape::*_changed() -> SunShape::changed()
    connect(this, &SunShape::shape_changed, this, &SunShape::changed);
    connect(this, &SunShape::sigma_changed, this, &SunShape::changed);
    connect(this, &SunShape::half_width_changed, this, &SunShape::changed);
    connect(this, &SunShape::csr_changed, this, &SunShape::changed);
    connect(m_custom_distribution,
            &SunShapeModel::changed,
            this,
            &SunShape::changed);

    // Initialization
    regenerate();
    update_current_distribution();
}

void SunShape::reset_current_distribution() {
    custom_distribution()->reset({
        { .angle = 0, .intensity = 1 },
        { .angle = 1, .intensity = 0.9 },
        { .angle = 2, .intensity = 0 },
    });
}

void SunShape::regenerate() {
    switch (m_shape) {
    case Shape::Gaussian: sample_gaussian(); break;
    case Shape::Pillbox: sample_pillbox(); break;
    case Shape::Buie_CSR: sample_buie(); break;
    case Shape::Custom: m_generated_distribution->clear(); break;
    case Shape::LimbDarkened: sample_limb_darkened(); break;
    }
    update_x_axis();
}

void SunShape::sample_gaussian() {
    // Point generation code referenced from app/src/sunshape
    // (SunShapeForm::UpdatePlot())

    int    num_points = 100;
    double theta_x    = 0;
    double theta_inc  = 3 * m_sigma / num_points;

    QVector<SunShapePoint> points;
    points.reserve(num_points);

    for (int i = 0; i < num_points; i++) {
        points.push_back({
            .angle     = theta_x,
            .intensity = 1.0 / exp(theta_x * theta_x / (2 * m_sigma * m_sigma)),
        });
        theta_x += theta_inc;
    }

    m_generated_distribution->reset(points);
}

void SunShape::sample_pillbox() {
    // Point generation code referenced from app/src/sunshape
    // (SunShapeForm::UpdatePlot())

    m_generated_distribution->reset({
        { .angle = 0, .intensity = 1 },
        { .angle = m_half_width, .intensity = 1 },
        { .angle = m_half_width + 1e-4, .intensity = 0 },
    });
}

void SunShape::sample_buie() {
    constexpr double min_csr_exclusive      = 0.0;
    constexpr double max_supported_csr      = 0.8;
    constexpr double solar_disk_radius_mrad = 4.65;
    constexpr double max_sample_angle_mrad  = 43.6;
    constexpr double sample_step_mrad       = 0.01;
    constexpr int    sample_count_estimate  = 4361;

    if (m_csr <= min_csr_exclusive || m_csr > max_supported_csr) {
        m_generated_distribution->reset();
        return;
    }

    // Buie's circumsolar model is driven by CSR, but the published intensity
    // equations use chi. These piecewise fits map the supported CSR range into
    // chi before deriving the aureole power-law terms below.
    const auto chi_from_csr = [](double csr) {
        if (csr > 0.145) {
            return -0.04419909985804843 +
                   csr * (1.401323894233574 +
                          csr * (-0.3639746714505299 +
                                 csr * (-0.9579768560161194 +
                                        1.1550475450828657 * csr)));
        }

        if (csr > 0.035) {
            return 0.022652077593662934 +
                   csr *
                       (0.5252380349996234 +
                        (2.5484334534423887 - 0.8763755326550412 * csr) * csr);
        }

        return 0.004733749294807862 +
               csr * (4.716738065192151 +
                      csr * (-463.506669149804 +
                             csr * (24745.88727411664 +
                                    csr * (-606122.7511711778 +
                                           5521693.445014727 * csr))));
    };

    const auto disk_intensity = [](double theta_mrad) {
        return std::cos(0.326 * theta_mrad) / std::cos(0.308 * theta_mrad);
    };

    const double chi   = chi_from_csr(m_csr);
    const double kappa = 0.9 * std::log(13.5 * chi) * std::pow(chi, -0.3);
    const double gamma = 2.2 * std::log(0.52 * chi) * std::pow(chi, 0.43) - 0.1;
    const double disk_edge_intensity = disk_intensity(solar_disk_radius_mrad);

    QVector<SunShapePoint> points;
    points.reserve(sample_count_estimate);

    // Store the non-negative radial profile. The graph mirrors it on demand
    // for display, while backend exports keep the radial profile unchanged.
    for (double theta = 0.0; theta <= max_sample_angle_mrad;
         theta += sample_step_mrad) {
        double intensity;

        if (theta <= solar_disk_radius_mrad) {
            intensity = disk_intensity(theta);
        } else {
            // Beyond the disk edge, Buie models the circumsolar aureole as a
            // power law. Cap it at the disk-edge value to avoid a
            // discontinuity.
            intensity = std::exp(kappa) * std::pow(theta, gamma);
            intensity = std::min(intensity, disk_edge_intensity);
        }

        points.push_back({
            .angle     = theta,
            .intensity = intensity,
        });
    }

    m_generated_distribution->reset(points);
}

void SunShape::sample_limb_darkened() {
    constexpr double disk_edge  = 4.65;
    constexpr int    num_points = 100;
    double           theta      = 0.0;
    double           theta_inc  = disk_edge / num_points;

    QVector<SunShapePoint> points;
    points.reserve(num_points + 1);

    for (int i = 0; i <= num_points; ++i) {
        points.push_back({
            .angle     = theta,
            .intensity = std::cos(0.326 * theta) / std::cos(0.308 * theta),
        });
        theta += theta_inc;
    }

    points.push_back({
        .angle = disk_edge + 1.e-3,
        .intensity = 0.0,
    });

    m_generated_distribution->reset(points);
}

void SunShape::update_x_axis() {
    QPointer<SunShapeModel> gdist = m_generated_distribution;
    QPointer<SunShapeModel> cdist = m_custom_distribution;

    switch (m_shape) {
    case Shape::Gaussian:
        gdist->set_x_axis_from(-3.3 * m_sigma);
        gdist->set_x_axis_to(3.3 * m_sigma);
        break;
    case Shape::Pillbox:
        gdist->set_x_axis_from(-3.3 * m_half_width);
        gdist->set_x_axis_to(3.3 * m_half_width);
        break;
    case Shape::Buie_CSR:
        gdist->set_x_axis_from(-20.0);
        gdist->set_x_axis_to(20.0);
        break;
    case Shape::LimbDarkened:
        gdist->set_x_axis_from(-1.3 * 4.65);
        gdist->set_x_axis_to(1.3 * 4.65);
        break;
    case Shape::Custom:
        // Code referenced from app/src/sunshape (SunShapeForm::UpdatePlot())
        if (cdist->count() >= 2) {
            double max_x = 0;
            for (int i = 0; i < cdist->count(); i++) {
                double angle = std::abs(cdist->get_at(i)->angle);
                if (angle > max_x) max_x = angle;
            }

            if (max_x == 0) {
                cdist->set_x_axis_from(-1.3);
                cdist->set_x_axis_to(1.3);
            } else {
                cdist->set_x_axis_from(-1.3 * max_x);
                cdist->set_x_axis_to(1.3 * max_x);
            }
        }
        break;
    }
}

void SunShape::update_current_distribution() {
    if (m_shape == Shape::Custom)
        set_current_distribution(m_custom_distribution);
    else
        set_current_distribution(m_generated_distribution);
}


Data::SunShape SunShape::get_sunshape_data() const {
    switch (m_shape) {
    case Shape::Gaussian: return Data::SunShape::GAUSSIAN;
    case Shape::Pillbox: return Data::SunShape::PILLBOX;
    case Shape::Buie_CSR: return Data::SunShape::BUIE_CSR;
    case Shape::Custom: return Data::SunShape::USER_DEFINED;
    case Shape::LimbDarkened: return Data::SunShape::LIMBDARKENED;
    default: return Data::SunShape::UNKNOWN;
    }
}


void SunShape::export_custom_distribution(QUrl url) {
    auto path = url.toLocalFile();

    auto file = QFile(path);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        emit notify(ANotification::error(
            QStringLiteral("Could not write CSV file. Check folder permissions "
                           "and available disk space.")));
        return;
    }

    auto stream = QTextStream(&file);

    stream << "angle,intensity\n";

    for (auto const& point : (*m_custom_distribution)) {
        stream << point.angle << "," << point.intensity << '\n';
    }
}

void SunShape::import_custom_distribution(QUrl url) {
    auto path = url.toLocalFile();

    auto file = QFile(path);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit notify(ANotification::error(
            QStringLiteral("Could not open CSV file. Check permissions.")));
        return;
    }

    auto stream = QTextStream(&file);

    auto first = stream.readLine();

    auto angle_idx     = first.indexOf("angle", 0, Qt::CaseInsensitive);
    auto intensity_idx = first.indexOf("intensity", 0, Qt::CaseInsensitive);

    auto angle_loc     = angle_idx > intensity_idx;
    auto intensity_loc = angle_idx < intensity_idx;

    QList<SunShapePoint> new_custom_shape;

    bool problems_detected = false;

    auto line_parser = [&](QStringView line) {
        auto parts = QStringView(line).split(',', Qt::SkipEmptyParts);

        bool ok = false;

        auto angle = parts.value(angle_loc).trimmed().toDouble(&ok);
        if (!ok) { problems_detected = true; }

        auto intensity = parts.value(intensity_loc).trimmed().toDouble(&ok);

        if (!ok) { problems_detected = true; }

        new_custom_shape << SunShapePoint {
            .angle     = parts.value(angle_loc).toDouble(),
            .intensity = parts.value(intensity_loc).toDouble()
        };
    };

    if (angle_idx < 0 or intensity_idx < 0) {
        angle_loc     = 0;
        intensity_loc = 1;

        line_parser(first.trimmed());
    }

    while (!stream.atEnd()) {
        auto line = stream.readLine();

        if (line.isEmpty()) continue;

        line_parser(line);
    }

    if (problems_detected) {
        emit notify(ANotification::error(
            QStringLiteral("Some values of the CSV could not be interpreted. "
                           "Check formatting.")));

        return;
    }

    m_custom_distribution->reset(new_custom_shape);
}

} // namespace SolTrace::GUI::App
