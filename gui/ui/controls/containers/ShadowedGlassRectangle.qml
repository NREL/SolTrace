import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects

Item {
    id: control

    property alias glassColor: background.glassColor

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

    GlassRectangle {
        id: background
        blur_source: control.blur_source

        anchors.fill: parent

        radius: control.radius
    }
}
