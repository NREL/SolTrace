import QtCore
import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts
import QtQuick.Dialogs

import SolTrace

ApplicationWindow {
    id: main_window

    title: "SolTrace - " + AppData.current_version_info

    width: 1280
    height: 720
    visible: true

    minimumWidth: 700
    minimumHeight: 600

    Material.theme: Material.Dark
    Material.accent: Material.Blue
    Material.foreground: App.theme.fontColor

    font.pointSize: App.theme.labelSize
    font.family: "Roboto"
    font.features: {"liga" : 0, "clig" : 0}

    SimulationScene {
        id: simulation_scene

        anchors.fill: parent
    }

    CorePanel {
        blur_source: simulation_scene
        anchors.fill: parent
    }
}
