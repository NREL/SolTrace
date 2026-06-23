import QtQuick
import QtQuick.Controls.Material
import SolTrace

ToolTip {
    id: root
    delay: 500
    contentItem: Label {
        text: root.text
        color: App.theme.fontColor
    }
}
