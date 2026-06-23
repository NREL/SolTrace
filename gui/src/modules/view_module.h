#pragma once

#include "utilities/qt_helpers.h"
#include <QColor>
#include <QObject>
#include <QQmlEngine>

namespace SolTrace::GUI::App {

class PanelData : public QObject {
    Q_OBJECT
    QML_ELEMENT

    const inline static QVector<int> m_sizes      = { 250, 550, 750, 9999 };
    const inline static QVector<int> m_thresholds = { 400, 600, 850 };

public:
    explicit PanelData(QObject* parent = nullptr);

    enum PanelSize { Small = 0, Normal = 1, Wide = 2, Full = 3 };

    Q_ENUM(PanelSize)

    Q_PROPERTY(QVector<int> sizes READ sizes CONSTANT)
    QVector<int> sizes() const;
    Q_PROPERTY(QVector<int> thresholds READ thresholds CONSTANT)
    QVector<int> thresholds() const;

    // Source of Truth
    Q_WRITABLE_PROPERTY(int, width, m_sizes[PanelSize::Normal])

    // Derive PanelSize enum from width
    Q_WRITABLE_PROPERTY(PanelSize, size, PanelSize::Normal)

    Q_WRITABLE_PROPERTY(bool, visible, false)
    Q_WRITABLE_PROPERTY(bool, saved_visible, false)

    Q_WRITABLE_PROPERTY(bool, inline_docs, false)

    // Q_WRITABLE_PROPERTY(bool, tags, false)
    // Idea: show walkthrough tags just for this section

public slots:
    bool is_small();
    bool is_normal();
    bool is_wide();
    void update_size();
};
class SimulationViewState : public QObject {
    Q_OBJECT
    QML_ELEMENT

public:
    enum class Camera { WASD = 0, Orbital = 1 };

    enum class Perspective { Normal = 0, Orthographic = 1 };

    Q_ENUM(Camera)
    Q_ENUM(Perspective)

    Q_WRITABLE_PROPERTY(Camera, camera, Camera::Orbital)
    Q_WRITABLE_PROPERTY(Perspective, perspective, Perspective::Normal)

    Q_WRITABLE_PROPERTY(bool, sun_viz, true)
    Q_WRITABLE_PROPERTY(bool, blueprint_mode, false)
    Q_WRITABLE_PROPERTY(double, sun_viz_scale, 50)
    Q_WRITABLE_PROPERTY(QColor, sun_color, "yellow")
    Q_WRITABLE_PROPERTY(QColor, geometry_color, "white")
};

class ViewModule : public QObject {
    Q_OBJECT
    QML_ELEMENT

    // Helper
    bool shrink_panel(const QVector<int>& sizes, QPointer<PanelData>& p);

public:
    explicit ViewModule(QObject* parent = nullptr);

    enum MouseMode {
        Camera = 0,
        SelectElement = 1,
        SelectMaterial = 2,
        SelectGeometry = 3,
        EditElement = 4,
        PickRay = 5
    };

    Q_ENUM(MouseMode)

    // Panel Data
    QOBJECT_READONLY_PROPERTY(PanelData, left_panel)
    QOBJECT_READONLY_PROPERTY(PanelData, right_panel)
    QOBJECT_READONLY_PROPERTY(PanelData, settings_panel)

    // Left Panel Section State
    Q_WRITABLE_PROPERTY(int, workflow_phase, 0)
    Q_WRITABLE_PROPERTY(int, configure_section, 0)
    Q_WRITABLE_PROPERTY(int, simulate_section, 0)
    Q_WRITABLE_PROPERTY(int, analyze_section, 0)

    Q_WRITABLE_PROPERTY(int, sun_section, 0)

    // Right Panel Section State
    Q_WRITABLE_PROPERTY(int, right_panel_section, 0)

    // Editor State
    Q_WRITABLE_PROPERTY(bool, editing_material, false)
    Q_WRITABLE_PROPERTY(bool, editing_geometry, false)
    Q_WRITABLE_PROPERTY(bool, editing_layout, false)
    Q_WRITABLE_PROPERTY(bool, editing_appearance, false)

    // Settings State
    Q_WRITABLE_PROPERTY(int, documentation_section, 0)

    Q_WRITABLE_PROPERTY(bool, simulation_content_view, false)
    Q_WRITABLE_PROPERTY(bool, show_intersections, true)
    Q_WRITABLE_PROPERTY(MouseMode, mouse_mode, MouseMode::Camera)

    // Viewport State
    QOBJECT_READONLY_PROPERTY(SimulationViewState, sim)

public slots:
    void fit_panels(int  available_width,
                    bool expanding_right_panel = false,
                    bool resizing_window       = false,
                    int  margin                = 30);
};

} // namespace SolTrace::GUI::App
