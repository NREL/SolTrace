#include "modules/sun/sun_shape_model.h"

#include <QClipboard>
#include <QGuiApplication>
#include <QRegularExpression>

#include <algorithm>
#include <cmath>

namespace SolTrace::GUI::App {

namespace {

QVector<SunShapePoint> normalized_radial_points(QVector<SunShapePoint> points) {
    for (auto& point : points) {
        point.angle = std::abs(point.angle);
    }

    std::sort(points.begin(), points.end(), [](auto const& a, auto const& b) {
        return a.angle < b.angle;
    });

    QVector<SunShapePoint> merged;
    for (auto const& point : points) {
        if (!merged.empty() &&
            std::abs(merged.back().angle - point.angle) < 1.0e-9) {
            merged.back().intensity =
                std::max(merged.back().intensity, point.intensity);
        } else {
            merged.push_back(point);
        }
    }

    return merged;
}

} // namespace

SunShapeModel::SunShapeModel(QObject* parent) : StructTableModel(parent) {
    connect(this, &SunShapeModel::dataChanged, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::rowsInserted, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::rowsRemoved, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::modelReset, this, &SunShapeModel::changed);
}

std::vector<double> SunShapeModel::get_angle_data() {
    // Backend sun-shape data is radial-only. Views that need a symmetric curve
    // should mirror these non-negative points at presentation time.
    auto points = normalized_radial_points(m_records);

    std::vector<double> result;
    for (auto const& point : std::as_const(points)) {
        result.push_back(point.angle);
    }
    return result;
}

std::vector<double> SunShapeModel::get_intensity_data() {
    // Keep this paired with get_angle_data(): both export the same normalized,
    // non-negative radial profile to the backend.
    auto points = normalized_radial_points(m_records);

    std::vector<double> result;
    for (auto const& point : std::as_const(points)) {
        result.push_back(point.intensity);
    }
    return result;
}

QVariantList SunShapeModel::variant_data() {
    auto points = normalized_radial_points(m_records);

    QVariantList custom_shape;
    for (auto const& source : std::as_const(points)) {
        QVariantMap point;
        point["angle"]     = source.angle;
        point["intensity"] = source.intensity;
        custom_shape.append(point);
    }
    return custom_shape;
}

void SunShapeModel::set_variant_data(QVariantList data) {
    QVector<SunShapePoint> points;
    for (const auto& item : data) {
        QVariantMap point = item.toMap();
        points.push_back({
            .angle     = std::abs(point["angle"].toDouble()),
            .intensity = point["intensity"].toDouble(),
        });
    }
    reset(normalized_radial_points(points));
}

int SunShapeModel::count() const {
    return rowCount();
}

void SunShapeModel::append(double angle, double intensity) {
    StructTableModel::append({ std::abs(angle), intensity });
    emit countChanged();
}

void SunShapeModel::reset(QVector<SunShapePoint> points) {
    StructTableModel<SunShapePoint>::reset(points);
    emit countChanged();
}

void SunShapeModel::remove(int index) {
    if (index < 0 || index >= m_records.count()) return;
    remove_at(index);
    emit countChanged();
}

void SunShapeModel::clear() {
    reset();
}

void SunShapeModel::copy_to_clipboard() {
    QString text   = "Angle (mrad)\tIntensity\n";
    auto    points = normalized_radial_points(m_records);

    for (auto const& point : std::as_const(points)) {
        text += QString::number(point.angle) + "\t" +
                QString::number(point.intensity) + "\n";
    }
    QGuiApplication::clipboard()->setText(text);
}

void SunShapeModel::paste_from_clipboard() {
    QVariantList rows;
    QString      text = QGuiApplication::clipboard()->text();
    for (auto const& line : text.split('\n')) {
        if (line.trimmed() == "") continue;
        QStringList v = line.split(QRegularExpression("[\\t,]"));
        if (v.length() >= 2) {
            QVariantMap row;
            row["angle"]     = v[0].trimmed().toDouble();
            row["intensity"] = v[1].trimmed().toDouble();
            rows.append(row);
        }
    }
    set_variant_data(rows);
}


} // namespace SolTrace::GUI::App
