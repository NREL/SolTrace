import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root

    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth
    clip: true

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Scene"
        }

        Label {
            text: "Edit scene visualization options"
            wrapMode: Text.WordWrap
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Geometry"
            collapsed: false

            ColorPickerField {
                id: geometryColorPicker
                color: App.view.sim.geometry_color
                label: "Global Geometry Color"
                onUpdated: {
                    App.view.sim.geometry_color = geometryColorPicker.color
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Sky"
            collapsed: false

            RowLayout {
                Layout.fillWidth: true


                STComboBox {
                    model: ["Day", "Blueprint", "Adaptive"]
                    currentIndex: App.view.sim.sky
                    onCurrentIndexChanged: App.view.sim.sky = currentIndex
                }

                STIconButton {
                    icon: "\uf059"
                    onClicked: skyInfo.open()

                    STPopup {
                        id: skyInfo

                        width: 300
                        y: -height - 8

                        Label {
                            width: parent.width
                            anchors.centerIn: parent
                            anchors.margins: 10
                            textFormat: Text.MarkdownText
                            wrapMode: Text.WordWrap
                            text: "**Set the environment sky to one of the following options:**\n" +
                                  "- **Adaptive**: adapts to ray source elevation \n" +
                                  "- **Blueprint**: neutral gray gradient\n" +
                                  "- **Day**: bright blue sky"
                        }

                    }
                }
            }
        }


        STPropertyPanel {
            Layout.fillWidth: true

            title: "Ray Source"
            collapsed: false

            ColumnLayout {
                Layout.fillWidth: true

                ColorPickerField {
                    id: sunColorPicker
                    label: "Color"
                    color: App.view.sim.sun_color
                    onUpdated: App.view.sim.sun_color = sunColorPicker.color
                }


                STSpinBoxField {
                    value: App.view.sim.sun_viz_scale
                    onValueModified: App.view.sim.sun_viz_scale = value
                    from: 0
                    to: 100
                    label: "Scale Factor"
                }

                STSwitch {
                    text: "Visibility"
                    checked: App.view.sim.sun_viz
                    onToggled: App.view.sim.sun_viz = checked
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
