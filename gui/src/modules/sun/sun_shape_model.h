#pragma once

#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

#include <QObject>
#include <QVariantList>
#include <QVector>

#include <vector>

namespace SolTrace::GUI::App {

/// One angular sample in a sun-shape intensity distribution.
struct SunShapePoint {
    double angle     = 0.0;
    double intensity = 0.0;

    RECORD_META(SunShapePoint, SM_EXPOSE_RW(angle), SM_EXPOSE_RW(intensity))
};

/// Editable table model for sun-shape samples.
class SunShapeModel : public StructTableModel<SunShapePoint> {
    Q_OBJECT
public:
    explicit SunShapeModel(QObject* parent = nullptr);

    /// Copy angle column values into a standard vector.
    std::vector<double> get_angle_data();

    /// Copy intensity column values into a standard vector.
    std::vector<double> get_intensity_data();

    Q_WRITABLE_PROPERTY(double, x_axis_from, -5)
    Q_WRITABLE_PROPERTY(double, x_axis_to, 5)
    Q_WRITABLE_PROPERTY(double, y_axis_from, 0)
    Q_WRITABLE_PROPERTY(double, y_axis_to, 1.2)

    /// Return all samples as QVariant values for QML serialization.
    QVariantList variant_data();

    /// Replace samples from QML-provided QVariant values.
    void set_variant_data(QVariantList data);

public slots:
    /// Replace all samples and update graph axis bounds.
    void reset(QVector<SunShapePoint> points = { });

    /// Append one sample point.
    void append(double angle = 0.0, double intensity = 0.0);

    /// Remove the sample at index.
    void remove(int index);

    /// Remove all samples.
    void clear();

    /// Copy samples to the system clipboard as text.
    void copy_to_clipboard();

    /// Paste samples from the system clipboard.
    void paste_from_clipboard();

    /// Number of sample rows.
    int count() const;

signals:
    void countChanged();
    void changed();
};

} // namespace SolTrace::GUI::App
