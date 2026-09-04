import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Rectangle {
    id: control

    signal clicked

    property alias icon: icon_label.text
    property string label: ""
    property real padding: 10
    property real iconSize: 16
    property real spacing: 6
    property string toolTip: ""
    readonly property alias containsMouse: mouse_area.containsMouse
    readonly property bool hasLabel: label.length > 0

    width: implicitWidth
    height: implicitHeight

    radius: height / 2

    color: {
        if (mouse_area.containsPress) return Material.highlightedRippleColor
        if (mouse_area.containsMouse) return Material.rippleColor
        return "transparent"
    }

    border.color: Material.dividerColor
    border.width: mouse_area.containsPress ? 1 : 0

    implicitHeight: Math.max(icon_label.implicitHeight,
                             text_label.implicitHeight) + padding
    implicitWidth: hasLabel ?
                       icon_label.implicitWidth + spacing
                       + text_label.implicitWidth + padding
                     : implicitHeight

    Row {
        anchors.centerIn: parent
        spacing: control.hasLabel ? control.spacing : 0

        Label {
            id: icon_label

            anchors.verticalCenter: parent.verticalCenter

            font.family: "Font Awesome 7 Free"
            font.pointSize: control.iconSize

            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter
        }

        Label {
            id: text_label

            anchors.verticalCenter: parent.verticalCenter

            visible: control.hasLabel
            text: control.label
            font.pointSize: App.theme.labelSize
            verticalAlignment: Qt.AlignVCenter
        }
    }

    MouseArea {
        id: mouse_area
        anchors.fill: control
        hoverEnabled: true
        onClicked: control.clicked()
    }

    STToolTip {
        visible: control.containsMouse && control.toolTip.length > 0
        text: control.toolTip
    }
}
