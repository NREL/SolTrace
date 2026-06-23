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
            text: "Theme"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Glass"
            collapsed: false

            ColumnLayout {
                ColorPickerField {
                    id: glassColorPicker
                    color: App.theme.glassColor
                    label: "Glass Color"
                    onUpdated: {
                        App.theme.glassColor.r = glassColorPicker.color.r
                        App.theme.glassColor.g = glassColorPicker.color.g
                        App.theme.glassColor.b = glassColorPicker.color.b
                    }
                }

                SliderField {
                    Layout.preferredWidth: 250
                    Layout.maximumWidth: 350

                    from: 0
                    to: 100
                    value: App.theme.glassColor.a * 100
                    onModified: {
                        App.theme.glassColor.a = value / 100
                    }

                    text: "Glass Translucency"
                }

                STDangerousButton {
                    Layout.preferredWidth: 250
                    text: "Reset"

                    onClicked: {
                        App.theme.glassColor.r = 0
                        App.theme.glassColor.g = 0
                        App.theme.glassColor.b = 0
                        App.theme.glassColor.a = 0.15
                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Font"
            collapsed: false

            ColumnLayout {
                ColorPickerField {
                    id: fontColorPicker

                    label: "Font Color"
                    onUpdated: {
                        App.theme.fontColor.r = fontColorPicker.color.r
                        App.theme.fontColor.g = fontColorPicker.color.g
                        App.theme.fontColor.b = fontColorPicker.color.b
                    }
                }

                STSpinBoxField {
                    label: "Header Size"
                    from: 1
                    to: 100
                    Layout.preferredWidth: 200
                    value: App.theme._headerSize
                    onValueModified: {
                        App.theme._headerSize = value
                    }
                }

                STSpinBoxField {
                    label: "Subheader Size"
                    from: 1
                    to: 100
                    Layout.preferredWidth: 200
                    value: App.theme._subHeaderSize
                    onValueModified: {
                        App.theme._subHeaderSize = value
                    }
                }

                STSpinBoxField {
                    label: "Label Size"
                    from: 1
                    to: 100
                    Layout.preferredWidth: 200
                    value: App.theme._labelSize
                    onValueModified: {
                        App.theme._labelSize = value
                    }
                }

                STSpinBoxField {
                    label: "Normal Font Size"
                    from: 1
                    to: 100
                    Layout.preferredWidth: 200
                    value: App.theme._normalSize
                    onValueModified: {
                        App.theme._normalSize = value
                    }
                }

                SliderField {
                    Layout.preferredWidth: 250
                    Layout.maximumWidth: 350

                    from: 0
                    to: 200
                    value: App.theme.zoomLevel * 100
                    onModified: {
                        App.theme.zoomLevel = value / 100
                    }

                    text: "Zoom Level"
                }

                STDangerousButton {
                    Layout.preferredWidth: 250
                    text: "Reset"

                    onClicked: {
                        App.theme.fontColor.r = 1
                        App.theme.fontColor.g = 1
                        App.theme.fontColor.b = 1
                        App.theme._headerSize = 17
                        App.theme._subHeaderSize = 16
                        App.theme._labelSize = 13
                        App.theme._normalSize = 15
                        App.theme.zoomLevel = 1
                    }
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
