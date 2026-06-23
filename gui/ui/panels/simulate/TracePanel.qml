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

    property bool singleColumn: App.view.left_panel.size === PanelData.Small
    property var labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    ColumnLayout {
        width: root.availableWidth

        InlineDocumentation {
            key: "placeholder_small"
            target: AppData.view.left_panel
            title: "Simulation Runner"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            collapsible: false
            title: "New Job"
            columns: root.singleColumn ? 1 : 2

            STComboBox {
                Layout.fillWidth: true
                Layout.columnSpan: root.singleColumn ? 1 : 2
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
        }
    }
}
