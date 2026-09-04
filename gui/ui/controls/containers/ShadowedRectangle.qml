import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects

Item {
    id: control

    property color glassColor: App.theme.glassColor

    property var blur_source

    property real radius: 10

    default property alias content: background.children

    RectangularShadow {
        anchors.fill: background
        offset.x: 0
        offset.y: 0
        radius: background.radius
        blur: 30
        spread: 3
        color: Material.dropShadowColor
    }

    Rectangle {
        id: background
        color: control.glassColor
        border.width: 1
        border.color: Qt.rgba(0.9, 0.9, 0.9, 0.15)
        anchors.fill: parent
        radius: control.radius
    }
}
