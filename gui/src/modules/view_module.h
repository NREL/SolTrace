#pragma once

#include "utilities/qt_helpers.h"
#include <QColor>
#include <QObject>
#include <QQmlEngine>

namespace SolTrace::GUI::App {

class FullPanelData : public QObject {
    Q_OBJECT
    QML_ELEMENT
public:
    explicit FullPanelData(QObject* parent = nullptr);

    enum class FullPanelMode {
        Settings         = 0,
        Documentation    = 1,
        BuildInformation = 2
    };

    Q_ENUM(FullPanelMode)

    Q_WRITABLE_PROPERTY(bool, visible, false)
    Q_WRITABLE_PROPERTY(bool, saved_visible, false)
    Q_WRITABLE_PROPERTY(FullPanelMode, mode, FullPanelMode::Settings)
    Q_WRITABLE_PROPERTY(int, settings_section, 0)
    Q_WRITABLE_PROPERTY(int, docs_section, 0)
    Q_WRITABLE_PROPERTY(int, build_section, 0)

public slots:
    void show();
    void hide();
};

class SplitPanelData : public QObject {
    Q_OBJECT
    QML_ELEMENT

    const inline static QVector<int> m_sizes      = { 450, 550, 750, 9999 };
    const inline static QVector<int> m_thresholds = { 550, 750, 850 };

public:
    explicit SplitPanelData(QObject* parent = nullptr);

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

    // Q_WRITABLE_PROPERTY(bool, tags, false)
    // Idea: show walkthrough tags just for this section

public slots:
    bool is_small();
    bool is_normal();
    bool is_wide();
    void update_size();

    void save_visibility();
    void restore_visibility();

    void show();
    void hide();
};

class SimulationViewState : public QObject {
    Q_OBJECT
    QML_ELEMENT

public:
    enum class Camera { WASD = 0, Orbital = 1 };

    enum class Perspective { Normal = 0, Orthographic = 1 };

    enum class SkyType { Day, Blueprint, Adaptive, Realistic };
    enum class RealisticSkyType { Clear, PartlyCloudy, LowSun };

    explicit SimulationViewState(QObject* parent = nullptr)
        : QObject(parent) { }

    Q_ENUM(Camera)
    Q_ENUM(Perspective)
    Q_ENUM(SkyType)
    Q_ENUM(RealisticSkyType)

    Q_WRITABLE_PROPERTY(Camera, camera, Camera::Orbital)
    Q_WRITABLE_PROPERTY(Perspective, perspective, Perspective::Normal)

    Q_WRITABLE_PROPERTY(bool, sun_viz, true)
    Q_WRITABLE_PROPERTY(SkyType, sky, SkyType::Adaptive)
    Q_WRITABLE_PROPERTY(RealisticSkyType,
                        realistic_sky,
                        RealisticSkyType::Clear)

    Q_WRITABLE_PROPERTY(double, sun_viz_scale, 50)
    Q_WRITABLE_PROPERTY(QColor, sun_color, "yellow")
    Q_WRITABLE_PROPERTY(QColor, geometry_color, "white")
    Q_WRITABLE_PROPERTY(bool, show_grid, true)
    Q_WRITABLE_PROPERTY(double, geometry_thickness, 0.05)
    Q_WRITABLE_PROPERTY(unsigned, geometry_subdivision_scale, 2)

    Q_WRITABLE_PROPERTY(int, fps_walk_speed, 10)
};

class ViewModule : public QObject {
    Q_OBJECT
    QML_ELEMENT

    // Helper
    bool shrink_panel(const QVector<int>& sizes, QPointer<SplitPanelData>& p);

public:
    explicit ViewModule(QObject* parent = nullptr);

    enum WorkflowPhase { Start, Load, Configure, Simulate, Analyze };

    enum MouseMode {
        Camera         = 0,
        SelectElement  = 1,
        SelectMaterial = 2,
        SelectGeometry = 3,
        EditElement    = 4,
        PickRay        = 5
    };

    Q_ENUM(MouseMode)
    Q_ENUM(WorkflowPhase)

    // Panel Data
    QOBJECT_READONLY_PROPERTY(SplitPanelData, left_panel)
    QOBJECT_READONLY_PROPERTY(SplitPanelData, right_panel)
    QOBJECT_READONLY_PROPERTY(FullPanelData, full_panel)

    // Left Panel Section State
    Q_WRITABLE_PROPERTY(WorkflowPhase, workflow_phase, Configure)
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

    // Simulation State
    Q_WRITABLE_PROPERTY(bool, simulation_content_view, false)
    Q_WRITABLE_PROPERTY(bool, show_intersections, true)
    Q_WRITABLE_PROPERTY(MouseMode, mouse_mode, MouseMode::Camera)
    Q_WRITABLE_PROPERTY(bool, inline_docs, false)

    // Viewport State
    QOBJECT_READONLY_PROPERTY(SimulationViewState, sim)

public slots:
    void fit_panels(int  available_width,
                    bool expanding_right_panel = false,
                    bool resizing_window       = false,
                    int  margin                = 30);

    void open_full_panel();
    void close_full_panel(int available_width);
    void toggle_full_panel(int available_width);

    void open_left_panel(int available_width);
    void close_left_panel();
    void toggle_left_panel(int available_width);

    void open_right_panel(int availalbe_width);
    void close_right_panel();
    void toggle_right_panel(int available_width);
};

} // namespace SolTrace::GUI::App
