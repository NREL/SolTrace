import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts
import QtGraphs

import SolTrace

Flickable {
    id: root
    property var left_panel_size: App.view.left_panel.size
    property var sun_module : App.sun
    readonly property real radToDeg: 180 / Math.PI
    readonly property real degToRad: Math.PI / 180

    function clamp(value, lo, hi) {
        return Math.max(lo, Math.min(hi, value))
    }

    function directionAzimuth() {
        const x = App.sun.position.x
        const y = App.sun.position.y
        var angle = Math.atan2(x, y) * radToDeg
        while (angle < 0) angle += 360
        while (angle >= 360) angle -= 360
        return angle
    }

    function directionElevation() {
        const x = App.sun.position.x
        const y = App.sun.position.y
        const z = App.sun.position.z
        const length = Math.sqrt(x * x + y * y + z * z)
        if (length <= 1.0e-12) return 0
        return Math.asin(clamp(z / length, -1, 1)) * radToDeg
    }

    function setDirectionAngles(azimuth, elevation) {
        const az = azimuth * degToRad
        const el = elevation * degToRad
        App.sun.position.x = Math.cos(el) * Math.sin(az)
        App.sun.position.y = Math.cos(el) * Math.cos(az)
        App.sun.position.z = Math.sin(el)
        App.sun.position.from_calculator = false
    }

    contentWidth: width
    contentHeight: content_column.implicitHeight
    clip: true
    boundsBehavior: Flickable.StopAtBounds

    ColumnLayout {
        id: content_column
        width: root.width

        SunPreview {
            id: sun_preview
            Layout.fillWidth: true
            Layout.preferredHeight: 148

        }

        InlineDocumentation {
            key: "configure.sun"
            target: App.view.left_panel

            Layout.margins: 8
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Visualization"
            collapsible: true

            ColumnLayout {
                width: parent.width

                STSwitch {
                    text: "Visible"
                    checked: App.view.sim.sun_viz
                    onToggled: App.view.sim.sun_viz = checked
                }

                ColorPickerField {
                    id: sunColorPicker
                    Layout.preferredWidth: 200
                    color: App.view.sim.sun_color
                    label: "Sun Color"
                    onUpdated: {
                        App.view.sim.sun_color = sunColorPicker.color
                    }
                }

                STSpinBoxField {
                    Layout.preferredWidth: 200

                    label: "Scale Factor"
                    value: App.view.sim.sun_viz_scale
                    onValueModified: {
                        App.view.sim.sun_viz_scale = value
                    }
                    from: 0
                    to: 100
                }
            }
        }


        STComboBar {
            id: bar
            currentIndex: App.view.sun_section
            onCurrentIndexChanged: App.view.sun_section = currentIndex
            Layout.fillWidth: true
            collapseLabels: App.view.left_panel.size === PanelData.Small
            model: ["Type & Position", "Shape"]
            iconModel: ["\uf0eb", "\uf1fe"]
        }

        StackLayout {
            currentIndex: App.view.sun_section

            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true

                RowLayout {
                    Layout.fillWidth: true

                    Label {
                        text: "Type:"
                    }

                    STComboBox {
                        Layout.fillWidth: true
                        currentIndex: App.sun.type
                        onCurrentIndexChanged: App.sun.type = currentIndex
                        model: ["Directional Sun", "Point Source"]
                    }
                }

                InlineDocumentation {
                    Layout.fillWidth: true
                    Layout.margins: 8

                    key: "configure.sun.type." + ["directional", "point_source"][App.sun.type]
                    target: App.view.left_panel
                }


                ColumnLayout {
                    Layout.fillWidth: true
                    visible: App.sun.type === SunModule.Directional

                    STPropertyPanel {
                        Layout.fillWidth: true

                        title: "Directional Angles"
                        collapsible: false

                        ColumnLayout {
                            width: parent.width

                            GridLayout {
                                Layout.fillWidth: true
                                columns: App.view.left_panel.size == PanelData.Small ? 1 : 2

                                STSpinBoxField {
                                    id: azimuthField
                                    Layout.fillWidth: true
                                    Layout.preferredWidth: 100
                                    label: "Azimuth"
                                    value: root.directionAzimuth()
                                    from: 0
                                    to: 360
                                    decimals: 3
                                    suffix: "deg"
                                    onValueModified: root.setDirectionAngles(
                                                         value,
                                                         root.directionElevation())
                                }

                                STSpinBoxField {
                                    id: elevationField
                                    Layout.fillWidth: true
                                    Layout.preferredWidth: 100
                                    label: "Elevation"
                                    value: root.directionElevation()
                                    from: -90
                                    to: 90
                                    decimals: 3
                                    suffix: "deg"
                                    onValueModified: root.setDirectionAngles(
                                                         root.directionAzimuth(),
                                                         value)
                                }
                            }

                            STButton {
                                Layout.alignment: Qt.AlignRight
                                text: "Set by Calculator"
                                onClicked: calculatorDialog.openForCurrent()
                            }
                        }
                    }

                    SunCalculatorDialog {
                        id: calculatorDialog
                    }
                }

                STPropertyPanel {
                    Layout.fillWidth: true

                    title: "Manual Position"
                    collapsible: false
                    visible: App.sun.type === SunModule.PointSource

                    GridLayout {
                        width: parent.width

                        STSpinBoxField {
                            Layout.row: 0
                            Layout.column: 0

                            Layout.fillWidth: true
                            Layout.preferredWidth: 100
                            value: App.sun.position.x
                            onValueModified: {
                                App.sun.position.x = value
                                App.sun.position.from_calculator = false
                            }
                            label: "X"
                            from: -1000000
                            to: 1000000
                        }

                        STSpinBoxField {
                            Layout.row: App.view.left_panel.size < 1 ? 1 : 0
                            Layout.column: App.view.left_panel.size < 1 ? 0 : 1

                            Layout.fillWidth: true
                            Layout.preferredWidth: 100
                            label: "Y"
                            value: App.sun.position.y
                            onValueModified: {
                                App.sun.position.y = value
                                App.sun.position.from_calculator = false
                            }
                            from: -1000000
                            to: 1000000
                        }

                        STSpinBoxField {
                            Layout.row: App.view.left_panel.size < 1 ? 2 : 0
                            Layout.column: App.view.left_panel.size < 1 ? 0 : 2

                            Layout.fillWidth: true
                            Layout.preferredWidth: 100
                            label: "Z"
                            value: App.sun.position.z
                            onValueModified: {
                                App.sun.position.z = value
                                App.sun.position.from_calculator = false
                            }
                            from: -1000000
                            to: 1000000
                        }
                    }
                }
            }

            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true

                RowLayout {
                    Layout.fillWidth: true
                    Label {
                        text: "Shape Type:"
                    }

                    STComboBox {
                        Layout.fillWidth: true
                        currentIndex: App.sun.shape.shape
                        onCurrentIndexChanged: App.sun.shape.shape = currentIndex
                        //collapseLabels: App.view.left_panel.size === PanelData.Small
                        model: ["Gaussian", "Pillbox", "CSR", "Custom", "Limb Darkened"]
                        //iconModel: ["\uf1fe", "\uf0c8", "\uf192", "\uf55b"]
                    }
                }

                InlineDocumentation {
                    Layout.fillWidth: true
                    Layout.margins: 8
                    key: "configure.sun.shape." + ["gaussian", "pillbox", "csr", "custom", "limb_darkened"][App.sun.shape.shape]
                    target: App.view.left_panel
                }

                GridLayout {
                    rows: 2
                    columns: 2
                    Layout.margins: 12

                    StackLayout {
                        currentIndex: App.sun.shape.shape
                        Layout.row: App.view.left_panel.size < 2 ? 1 : 0
                        Layout.column: App.view.left_panel.size < 2 ? 0 : 1
                        Layout.alignment: Qt.AlignHCenter
                        Layout.preferredWidth: 200

                        STSpinBoxField {
                            Layout.fillWidth: true
                            label: "Standard Deviation"
                            value: App.sun.shape.sigma
                            decimals: 3
                            suffix: "mrad"
                            onValueChanged: { App.sun.shape.sigma = value }
                        }

                        STSpinBoxField {
                            Layout.fillWidth: true
                            label: "Half-width"
                            value: App.sun.shape.half_width
                            decimals: 3
                            suffix: "mrad"
                            onValueChanged: { App.sun.shape.half_width = value }
                        }

                        STSpinBoxField {
                            Layout.fillWidth: true
                            label: "Circumsolar Ratio"
                            value: App.sun.shape.csr
                            decimals: 3
                            from: 0
                            to: 0.8
                            stepSize: 0.1
                            onValueChanged: { App.sun.shape.csr = value }
                        }

                        CustomSunShapeTable { }

                        Label {
                            Layout.fillWidth: true
                            text: "No parameters"
                            horizontalAlignment: Text.AlignHCenter
                            opacity: 0.7
                        }
                    }

                    SunShapeGraph {
                        Layout.row: 0
                        Layout.column: 0
                        Layout.preferredWidth: 300
                        Layout.preferredHeight: 400
                        Layout.fillWidth: true
                        Layout.fillHeight: true

                        title: "Emission Profile"
                        xAxisTitle: "Angle (mrad)"
                        yAxisTitle: "Intensity"
                    }
                }
            }
        }
    }

    ScrollBar.vertical: ScrollBar {
        id: vbar
    }
}
