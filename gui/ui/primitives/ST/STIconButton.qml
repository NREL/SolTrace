import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Rectangle {
    id: control

    signal clicked

    property alias text: label.text
    property alias label: label
    property real padding: 10
    property real iconSize: 16
    property string toolTip: ""
    readonly property alias containsMouse: mouse_area.containsMouse

    width: label.width + padding
    height: width

    radius: width / 2

    color: mouse_area.containsPress ?
               Material.highlightedRippleColor : "transparent"

    border.color: Material.dividerColor
    border.width: mouse_area.containsPress ? 1 : 0

    implicitHeight: Math.max(label.implicitHeight, label.implicitWidth) + padding
    implicitWidth: implicitHeight


    Label {
        id: label
        font.family: "Font Awesome 7 Free"

        font.pointSize: control.iconSize

        
        anchors.centerIn: parent
    }

    // Rectangle {
    //     anchors.fill: parent
    //     radius: width / 2
    // }

    MouseArea {
        id: mouse_area
        anchors.fill: parent
        hoverEnabled: true
        onClicked: control.clicked()
    }

    STToolTip {
        visible: control.containsMouse && control.toolTip.length > 0
        text: control.toolTip
    }
}

