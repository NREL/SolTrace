import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Button {
    id: control

    padding: 0

    property string text_icon

    property color down_color: Material.rippleColor
    property color idle_color: App.theme.glassColor

    contentItem: RowLayout {
        opacity: enabled ? 1.0 : 0.3

        anchors.fill: parent
        spacing: 6

        Item {
            Layout.fillWidth: true
            Layout.fillHeight: true
            //color: "blue"
        }

        Label {
            visible: control.text_icon.length > 0
            text: control.text_icon
            horizontalAlignment: Text.AlignHCenter
            verticalAlignment: Text.AlignVCenter

            font.family: "Font Awesome 7 Free"

            //Layout.fillWidth: true
            Layout.fillHeight: true
        }

        Label {
            text: control.text
            horizontalAlignment: Text.AlignHCenter
            verticalAlignment: Text.AlignVCenter
            elide: Text.ElideRight


            

            //Layout.fillWidth: true
            Layout.fillHeight: true
        }

        Item {
            Layout.fillWidth: true
            Layout.fillHeight: true
            //color: "blue"
        }
    }

    background: Item {
        implicitWidth: 100
        implicitHeight: 36

        opacity: enabled ? 1 : 0.3

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
            anchors.fill: parent
            radius: height / 2
            color: control.down ?
                       control.down_color
                     :
                       control.idle_color


            border.color: Material.dividerColor

            border.width: 1
        }
    }
}
