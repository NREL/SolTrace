import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STPopup {

    margins: 10

    width: 320

    // Navigation
    ColumnLayout {
        anchors.fill: parent
        InlineDocumentation {
            key: "placeholder_small"
            title: "Cameras"
            target: AppData.view.left_panel
        }

        STComboBar {
            currentIndex: AppData.view.sim.camera
            onCurrentIndexChanged: AppData.view.sim.camera = currentIndex
            Layout.fillWidth: true
            model: ["FPS Camera", "Orbital Camera"]
            iconModel: ["\uf03d", "\uf135"]
        }

        InlineDocumentation {
            key: "placeholder_small"
            title: "Camera Perspectives"
            target: AppData.view.left_panel
        }

        STComboBar {
            currentIndex: AppData.view.sim.perspective
            onCurrentIndexChanged: AppData.view.sim.perspective = currentIndex
            Layout.fillWidth: true
            model: ["Perspective", "Orthographic"]
            iconModel: ["\uf1b2", "\uf0c8"]
        }

        GridLayout {
            columns: 2

            Layout.fillWidth: true

            uniformCellWidths: true
            uniformCellHeights: true

            Label {
                text: "Align View"
                Layout.columnSpan: 2
            }

            STButton {
                text: "Reset"
                Layout.columnSpan: 2
                Layout.fillWidth: true
                onClicked: simulation_scene.reset_camera_view()
            }

            STButton {
                text: "-X"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.X, true)
            }

            STButton {
                text: "X"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.X, false)
            }

            STButton {
                text: "-Y"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.Y, true)
            }

            STButton {
                text: "Y"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.Y, false)
            }

            STButton {
                text: "-Z"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.Z, true)
            }

            STButton {
                text: "Z"
                onClicked: simulation_scene.align_to_axis(CameraController.Axis.Z, false)
            }
        }

        STSwitch {
            Layout.fillWidth: true

            text: "Blueprint Mode"
            checked: App.view.sim.blueprint_mode

            onToggled: {
                App.view.sim.blueprint_mode = !App.view.sim.blueprint_mode
            }
        }
    }

}
