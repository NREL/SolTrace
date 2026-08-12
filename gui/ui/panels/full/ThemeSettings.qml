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

    function resetGlass() {
        App.theme.glassColor = App.theme.defaultGlassColor
        App.theme.fontColor = App.theme.defaultFontColor
        glassColorPicker.color = App.theme.glassColor
        glassAlphaSlider.value = App.theme.glassColor.a * 100
        fontColorPicker.color = App.theme.fontColor
    }

    function resetFont() {
        App.theme.fontColor = App.theme.defaultFontColor
        App.theme._headerSize = App.theme.defaultHeaderSize
        App.theme._subHeaderSize = App.theme.defaultSubHeaderSize
        App.theme._labelSize = App.theme.defaultLabelSize
        App.theme._normalSize = App.theme.defaultNormalSize
        App.theme.zoomLevel = 1

        fontColorPicker.color = App.theme.fontColor
        headerSizeField.value = App.theme._headerSize
        subHeaderSizeField.value = App.theme._subHeaderSize
        labelSizeField.value = App.theme._labelSize
        normalSizeField.value = App.theme._normalSize
        zoomLevelField.value = App.theme.zoomLevel * 100
    }

    function reset() {
        resetGlass()
        resetFont()
    }

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Theme"
        }

        STDangerousButton {
            Layout.preferredWidth: 100

            text: "Reset"
            down_color: App.theme.defaultGlassColor

            onClicked: root.reset()
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
                    id: glassAlphaSlider

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
                    Layout.preferredWidth: 100
                    text: "Reset"

                    down_color: App.theme.defaultGlassColor

                    onClicked: root.resetGlass()
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
                    color: App.theme.fontColor
                    label: "Font Color"
                    onUpdated: {
                        App.theme.fontColor.r = fontColorPicker.color.r
                        App.theme.fontColor.g = fontColorPicker.color.g
                        App.theme.fontColor.b = fontColorPicker.color.b
                    }
                }

                STSpinBoxField {
                    id: headerSizeField
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
                    id: subHeaderSizeField
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
                    id: labelSizeField
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
                    id: normalSizeField
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
                    id: zoomLevelField
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
                    Layout.preferredWidth: 100
                    text: "Reset"

                    down_color: App.theme.defaultGlassColor
                    onClicked: root.resetFont()
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
