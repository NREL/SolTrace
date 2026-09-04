import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ColumnLayout {
    id: root
    property int from: 0
    property int to: 100
    property int value: 50
    property string text: "Label"

    signal modified()

    width: 300
    height: slider.implicitHeight + label.implicitHeight
    spacing: 0

    Slider {
        id: slider

        Layout.fillWidth: true
        palette.accent: App.theme.fontColor

        from: root.from
        to: root.to
        value: root.value
        onValueChanged: {
            root.value = value
            root.modified()
        }
    }

    Label {
        id: label
        Layout.alignment: Qt.AlignHCenter
        text: root.text + ": " + root.value
    }
}
