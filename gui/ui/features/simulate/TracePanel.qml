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

    property bool showProgressStatus: false
    property bool simulationRunning: AppData.simulation.is_running

    onSimulationRunningChanged: {
        if (simulationRunning) {
            hideIdleStatusTimer.stop()
            showProgressStatus = true
        } else if (showProgressStatus) {
            hideIdleStatusTimer.restart()
        }
    }

    Timer {
        id: hideIdleStatusTimer
        interval: 3000
        repeat: false
        onTriggered: root.showProgressStatus = AppData.simulation.is_running
    }

    Component.onCompleted: showProgressStatus = simulationRunning

    ColumnLayout {
        width: root.availableWidth

        InlineDocumentation {
            key: "simulate.trace"
            title: "Simulation Runner"
        }

        STFormPanel {
            title: "Trace Configuration"

            STComboBox {
                Layout.fillWidth: true
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

            STFormHelp {
                key: "simulate.trace.rays"
            }

            STFormRow {
                label: "# of Rays"

                STSpinBox {
                    id: rayCountField
                    Layout.fillWidth: true
                    from: 1
                    to: AppData.simulation.max_ray_count
                    onValueModified: AppData.simulation.ray_count = value

                    Binding {
                        target: rayCountField
                        property: "value"
                        value: AppData.simulation.ray_count
                        restoreMode: Binding.RestoreBinding
                    }
                }
            }

            STFormRow {
                label: "Max # Rays Traced"

                STSpinBox {
                    Layout.fillWidth: true
                    from: 1
                    value: AppData.simulation.max_ray_count
                    to: 1000000000
                    onValueModified: AppData.simulation.update_max_ray_count(value)
                }
            }

            STFormRow {
                label: "# of CPU Cores"
                visible: AppData.simulation.runner < 2

                STSpinBox {
                    Layout.fillWidth: true
                    from: 1
                    value: AppData.simulation.max_threads
                    to: 64
                    onValueModified: AppData.simulation.max_threads = value
                }
            }

            STFormRow {
                label: "Seed Value"

                STSpinBox {
                    Layout.fillWidth: true
                    from: 1
                    value: AppData.simulation.seed_value
                    to: 10000000
                    onValueModified: AppData.simulation.seed_value = value
                }
            }

            STFormHelp {
                key: "simulate.trace.options"
            }

            STFormRow {
                label: "Options"

                ColumnLayout {
                    Layout.fillWidth: true

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
            }
        }

        STFormPanel {
            collapsible: true
            title: "Execution"

            STFormRow {
                label: "Progress"
                visible: root.showProgressStatus

                ProgressBar {
                    Layout.fillWidth: true
                    from: 0
                    to: 100
                    value: AppData.simulation.progress

                    enabled: AppData.simulation.is_running
                }
            }

            STFormRow {
                label: "Stage"
                visible: root.showProgressStatus

                Label {
                    Layout.fillWidth: true
                    text: AppData.simulation.is_running ?
                              AppData.simulation.current_stage : "Idle"
                }
            }

            STButton {
                Layout.fillWidth: true
                text: "Start Trace"
                left_text_icon: "\uf0da"
                onClicked: {
                    AppData.simulation.run()
                }
            }

            STFormHelp {
                key: "simulate.trace.results"
            }

            STIconButton {
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
