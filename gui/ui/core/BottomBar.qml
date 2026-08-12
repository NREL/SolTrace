import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Controls.Material

import SolTrace

RowLayout {
    id: bottom_bar

    required property int available_width
    readonly property int bar_height: 42
    required property var blur_source
    property bool collapsed: false
    property int normal_width: left_bottom_bar.width + workflow_bar.width + right_bottom_bar.width + spacing * 2

    height: workflow_bar.is_open ? workflow_large_pane.implicitHeight : bar_height
    anchors.left: parent.left
    anchors.right: parent.right
    anchors.bottom: parent.bottom
    anchors.margins: 10

    ShadowedGlassRectangle {
        id: left_bottom_bar
        Layout.preferredWidth: docGroup.implicitWidth + 40
        Layout.preferredHeight: bottom_bar.bar_height
        Layout.alignment: Qt.AlignBottom

        blur_source: bottom_bar.blur_source
        radius: height / 2
        glassColor: App.theme.glassColor

        RowLayout {
            id: docGroup
            height: parent.height
            anchors.horizontalCenter: parent.horizontalCenter

            STIconButton {
                icon: "\uf02d"
                toolTip: "Docs"

                onClicked: {
                    App.view.full_panel.mode = FullPanelData.Documentation
                    if (!App.view.full_panel.visible) App.view.toggle_full_panel(bottom_bar.width)
                }
            }

            STIconButton {
                icon: "\uf05a"
                toolTip: "Build Info"

                onClicked: {
                    App.view.full_panel.mode = FullPanelData.BuildInformation
                    if (!App.view.full_panel.visible) App.view.toggle_full_panel(bottom_bar.width)
                }
            }

        }
    }

    Item {
        Layout.fillWidth: true
    }

    WorkflowPaneSmall {
        id: workflow_bar_small

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: 42

        blur_source: bottom_bar.blur_source
        radius: bottom_bar.bar_height / 2
        glassColor: App.theme.glassColor

        visible: bottom_bar.collapsed
    }

    ShadowedGlassRectangle {
        id: workflow_bar

        property bool is_open: false

        Layout.preferredWidth: workflow_bar.is_open ? workflow_large_pane.implicitWidth : workflow_data_pane.implicitWidth
        Layout.preferredHeight: workflow_bar.is_open ? workflow_large_pane.implicitHeight : 42
        Layout.alignment: Qt.AlignBottom

        blur_source: bottom_bar.blur_source
        radius: bottom_bar.bar_height / 2
        glassColor: App.theme.glassColor

        visible: !bottom_bar.collapsed

        onIs_openChanged: {
            if (is_open) {
                workflow_data_pane.visible = false
                workflow_data_pane.opacity = 0
                workflow_large_pane.opacity = 0
                workflow_large_pane.visible = true
                workflow_large_fade.restart()
            } else {
                workflow_large_pane.visible = false
                workflow_large_pane.opacity = 0
                workflow_data_pane.opacity = 0
                workflow_data_pane.visible = true
                workflow_small_fade.restart()
            }
        }

        Item {
            NumberAnimation {
                id: workflow_small_fade
                target: workflow_data_pane
                property: "opacity"
                to: 1
                duration: 150
            }

            NumberAnimation {
                id: workflow_large_fade
                target: workflow_large_pane
                property: "opacity"
                to: 1
                duration: 150
            }
        }

        WorkflowPane {
            id: workflow_data_pane

            anchors.fill: parent

            opacity: 1
            visible: true

            onOpenClicked: workflow_bar.is_open = true
            blur_source: bottom_bar.blur_source
        }

        WorkflowPaneLarge {
            id: workflow_large_pane

            anchors.fill: parent

            opacity: 0
            visible: false

            onCloseClicked: workflow_bar.is_open = false
        }
    }

    Item {
        Layout.fillWidth: true
    }

    ShadowedGlassRectangle {
        id: right_bottom_bar
        Layout.preferredWidth: cameraGroup.implicitWidth + 40
        Layout.preferredHeight: bottom_bar.bar_height
        Layout.alignment: Qt.AlignBottom

        blur_source: bottom_bar.blur_source
        radius: height / 2
        glassColor: App.theme.glassColor

        RowLayout {
            id: cameraGroup
            height: parent.height
            anchors.horizontalCenter: parent.horizontalCenter

            STIconButton {
                icon: AppData.view.sim.camera === SimulationViewState.WASD ? "\uf047" : "\ue4bb"
                toolTip: "Camera Mode"

                onClicked: camera_popup.open()

                STPopup {
                    id: camera_popup

                    width: 250
                    y: -height - 16
                    x: (parent.width - width) / 2

                    ColumnLayout {
                        width: parent.width

                        InlineDocumentation {
                            key: "view.camera"
                            title: "Cameras"
                        }

                        STComboBar {
                            currentIndex: AppData.view.sim.camera
                            onCurrentIndexChanged: AppData.view.sim.camera = currentIndex
                            Layout.fillWidth: true
                            model: ["FPS Camera", "Orbital Camera"]
                            iconModel: ["\uf03d", "\uf135"]
                        }

                        SliderField {
                            value: AppData.view.sim.fps_walk_speed
                            from: 10
                            to: 100
                            text: "FPS Camera Speed"
                            visible: AppData.view.sim.camera == 0

                            onValueChanged: {
                                App.view.sim.fps_walk_speed = value
                            }
                        }

                        InlineDocumentation {
                            key: "view.perspective"
                            title: "Camera Perspectives"
                        }

                        STComboBar {
                            currentIndex: AppData.view.sim.perspective
                            onCurrentIndexChanged: AppData.view.sim.perspective = currentIndex
                            Layout.fillWidth: true
                            model: ["Perspective", "Orthographic"]
                            iconModel: ["\uf1b2", "\uf0c8"]
                        }
                    }
                }
            }

            STIconButton {
                icon: "\uf279"
                toolTip: "View Tools"

                onClicked: view_tools_popup.open()

                STPopup {
                    id: view_tools_popup

                    width: 250
                    y: -height - 16
                    x: (parent.width - width) / 2

                    ColumnLayout {
                        width: parent.width

                        STButton {
                            Layout.fillWidth: true
                            text: "Fit All in View"
                            left_text_icon: "\uf279"
                            onClicked: {
                                simulation_scene.fit_all_in_view()
                                view_tools_popup.close()
                            }
                        }
                    }
                }
            }

            STIconButton {
                icon: "\uf568"
                toolTip: "Alignment"

                onClicked: camera_alignment_popup.open()

                STPopup {
                    id: camera_alignment_popup
                    width: 250
                    y: -height - 16
                    x: (parent.width - width) / 2

                    GridLayout {
                        columns: 2

                        Layout.fillWidth: true

                        uniformCellWidths: true
                        uniformCellHeights: true

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
                }
            }

            STIconButton {
                icon: "\uf625"
                toolTip: "Performance"

                onClicked: performance_popup.open()

                STPopup {
                    id: performance_popup
                    width: 260
                    y: -height - 16
                    x: (parent.width - width) / 2

                    GridLayout {
                        width: parent.width
                        columns: 2

                        STPropertySeparator {
                            title: "Performance"
                        }

                        STPropertyLabel {
                            text: "Surface Width"
                        }

                        STDoubleSpinBox {
                            Layout.fillWidth: true
                            from: 0
                            to: 10
                            stepSize: 0.01
                            decimals: 2
                            value: App.view.sim.geometry_thickness
                            onValueModified: App.view.sim.geometry_thickness = value
                        }

                        STPropertyLabel {
                            text: "Subdivision Scale"
                        }

                        STSpinBox {
                            Layout.fillWidth: true
                            from: 1
                            to: 10
                            value: App.view.sim.geometry_subdivision_scale
                            onValueModified: App.view.sim.geometry_subdivision_scale = value
                        }

                        STSwitch {
                            Layout.fillWidth: true
                            Layout.columnSpan: 2
                            text: "Show Grid"
                            checked: App.view.sim.show_grid
                            onToggled: App.view.sim.show_grid = checked
                        }
                    }
                }
            }
        }
    }
}
