#pragma once

#include "modules/sun/sun_shape_model.h"
#include "utilities/notification.h"

#include "sun.hpp"

#include <QList>
#include <QObject>
#include <QUrl>
#include <qqmlintegration.h>

namespace SolTrace::GUI::App {

/// QML-facing sun-shape editor and sampler.
class SunShape : public QObject {
    Q_OBJECT
    QML_ELEMENT

    void regenerate();
    void sample_gaussian();
    void sample_pillbox();
    void sample_buie();
    void sample_limb_darkened();
    void update_x_axis();

    void update_current_distribution();

public:
    explicit SunShape(QObject* parent = nullptr);

    /// Convert current GUI settings to the library sun-shape object.
    SolTrace::Data::SunShape get_sunshape_data() const;

    // Note that this is separate from SolTrace's SunShape enum to maintain
    // independence from backend modifications
    enum class Shape { Gaussian, Pillbox, Buie_CSR, Custom, LimbDarkened };
    Q_ENUM(Shape)
    Q_WRITABLE_PROPERTY(Shape, shape, Shape::Gaussian)

    // Gaussian
    Q_WRITABLE_PROPERTY(double, sigma, 4.65)

    // Pillbox
    Q_WRITABLE_PROPERTY(double, half_width, 4.65)

    // Buie
    Q_WRITABLE_PROPERTY(double, csr, 0.1)

    // Generated distribution (using sigma, half_width, csr)
    QOBJECT_READONLY_PROPERTY(SunShapeModel, generated_distribution)

    // Custom distribution (using user-defined points)
    QOBJECT_READONLY_PROPERTY(SunShapeModel, custom_distribution)

    // Current distribution
    QOBJECT_WRITABLE_PROPERTY(SunShapeModel, current_distribution)

    /// Reset the active distribution to the generated/custom source for shape.
    void reset_current_distribution();

public slots:
    void export_custom_distribution(QUrl);
    void import_custom_distribution(QUrl);

signals:
    void changed();
    void notify(ANotification);
};

} // namespace SolTrace::GUI::App
