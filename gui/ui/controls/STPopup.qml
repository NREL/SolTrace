import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts
import QtQuick.Effects

Popup {
    id: root

    property bool enable_shadow: false

    background: Rectangle {
        radius: 14

        color: Qt.alpha(Material.backgroundColor, .80)

        RectangularShadow {
            enabled: root.enable_shadow
            anchors.fill: parent
            offset.x: 0
            offset.y: 0
            radius: parent.radius
            blur: 30
            spread: 2
            color: Material.dropShadowColor
        }
    }

    Overlay.modal: Rectangle {
        color: Qt.rgba(0, 0, 0, 0.5)
    }

    margins: 1
    padding: 8
}
