import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root

    Layout.fillHeight: true
    Layout.fillWidth: true

    contentWidth: availableWidth

    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property var labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    property int columnSpan: singleColumn ? 1 : 2

    ColumnLayout {
        width: root.availableWidth

        InlineDocumentation {
            key: "placeholder_small"
            target: AppData.view.left_panel
            title: "Simulation Runner"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Trace Configuration"
            columns: root.singleColumn ? 1 : 2

            STComboBox {
                Layout.fillWidth: true
                Layout.columnSpan: root.columnSpan
                model: AppData.simulation.runners
                textRole: "name"
                valueRole: "runner"
                currentIndex: AppData.simulation.runners.index_of(
                                  AppData.simulation.runner)
                onActivated: (index) => {
                    AppData.simulation.runner =
                            AppData.simulation.runners.runner_at(index)
                }
            }

            STPropertyLabel {
                text: "# of Rays"
                Layout.alignment: root.labelAlignment
            }

            STSpinBox {
                Layout.fillWidth: true
                from: 1
                value: AppData.simulation.ray_count
                to: 1000000000
                onValueModified: AppData.simulation.ray_count = value
            }

            STPropertyLabel {
                text: "Max # Rays Traced"
                Layout.alignment: root.labelAlignment
            }

            STSpinBox {
                Layout.fillWidth: true
                from: 1
                value: AppData.simulation.max_ray_count
                to: 1000000000
                onValueModified: AppData.simulation.max_ray_count = value
            }

            STPropertyLabel {
                text: "# of CPU Cores"
                Layout.alignment: root.labelAlignment
                visible: AppData.simulation.runner < 2
            }

            STSpinBox {
                Layout.fillWidth: true
                from: 1
                value: AppData.simulation.max_threads
                to: 64
                visible: AppData.simulation.runner < 2
                onValueModified: AppData.simulation.max_threads = value
            }

            STPropertyLabel {
                text: "Seed Value"
                Layout.alignment: root.labelAlignment
            }

            STSpinBox {
                Layout.fillWidth: true
                from: 1
                value: AppData.simulation.seed_value
                to: 10000000
                onValueModified: AppData.simulation.seed_value = value
            }

            STPropertyLabel {
                text: "Options"
                //Layout.columnSpan: root.columnSpan
                Layout.rowSpan: root.singleColumn ? 1 : 3
                Layout.alignment: root.labelAlignment
            }

            STSwitch {
                text: "Sun Shape"
                checked: AppData.simulation.sun_shape
                onToggled: AppData.simulation.sun_shape = checked
            }

            STSwitch {
                text: "Optical Errors"
                checked: AppData.simulation.optical_errors
                onToggled: AppData.simulation.optical_errors = checked
            }

            STSwitch {
                text: "Tracer Debug (if applicable)"
                //checked: AppData.simulation.optical_errors
                //onToggled: AppData.simulation.optical_errors = checked
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            collapsible: true
            title: "Execution"

            STPropertyLabel {
                text: "Progress"
            }

            ProgressBar {
                Layout.fillWidth: true
                from: 0
                to: 100
                value: AppData.simulation.progress

                enabled: AppData.simulation.is_running
            }

            STPropertyLabel {
                text: "Stage"
            }

            Label {
                Layout.fillWidth: true
                text: AppData.simulation.is_running ?
                          AppData.simulation.current_stage : "Idle"
            }

            STButton {
                Layout.columnSpan: 2
                Layout.fillWidth: true
                text: "Start Trace"
                left_text_icon: "\uf0da"
                onClicked: {
                    AppData.simulation.run()
                }
            }

            STIconButton {
                Layout.columnSpan: root.columnSpan
                Layout.fillWidth: false
                Layout.alignment: Qt.AlignRight
                icon: "\uf1da"
                label: "View Results"
                onClicked: {
                    App.view.workflow_phase = ViewModule.Analyze
                    App.view.analyze_section = 0
                    App.view.simulation_content_view = true
                }
            }
        }
    }
}
