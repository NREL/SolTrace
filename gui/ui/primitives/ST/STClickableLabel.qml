import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Label {
    id: control

    text: ""
    
    signal clicked

    property alias containsMouse: mouse_area.containsMouse

    property alias borderWidth: hover_rect.border.width

    Rectangle {
        id: hover_rect
        anchors.fill: parent
        radius: width / 2

        anchors.margins: -5

        color: mouse_area.containsMouse ? Material.rippleColor :
                   mouse_area.containsPress ? Material.highlightedRippleColor :
                                              "transparent"

        border.color: Material.dividerColor
        border.width: 0

        Behavior on color {
            ColorAnimation {
                duration: 100
            }
        }
    }

    MouseArea {
        id: mouse_area
        anchors.fill: parent
        onClicked: control.clicked()
        hoverEnabled: true
    }
}
