import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Switch {
    id: root

    property color down_color: Material.rippleColor
    property color idle_color: App.theme.glassColor

    padding: 2

    indicator: Rectangle {
        implicitWidth: 48
        implicitHeight: 26
        x: root.leftPadding
        y: parent.height / 2 - height / 2
        radius: 13
        color: root.checked ? Material.color(Material.Green) : root.idle_color
        border.color: Material.dividerColor

        RectangularShadow {
            anchors.fill: parent
            offset.x: 0
            offset.y: 0
            radius: parent.radius
            blur: 30
            spread: 3
            color: Material.dropShadowColor
        }

        Rectangle {
            x: root.checked ? parent.width - width : 0
            width: parent.implicitHeight
            height: width
            radius: 13
            color: root.down ? root.down_color : root.idle_color
            border.color: Material.rippleColor

            Behavior on x {
                NumberAnimation {
                    duration: 100
                }
            }

            Behavior on color {
                ColorAnimation {
                    duration: 100
                }
            }
        }
    }
}
