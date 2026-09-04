import QtQuick
import QtQuick.Controls.Material
import QtQuick.Dialogs
import QtQuick.Effects
import QtQuick.Layouts
import SolTrace

Item {
    id: root

    property color color: "white"
    property string label: ""

    width: 200
    height: column.implicitHeight

    signal updated()

    ColumnLayout {
        id: column

        width: parent.width

        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: 32

            RectangularShadow {
                anchors.fill: parent

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
                color: root.color

                border.color:  Material.dividerColor

                border.width: 1
            }
        }

        Label {
            id: label
            text: root.label
            Layout.fillWidth: true
            horizontalAlignment: Text.AlignHCenter
            visible: text.length > 0
        }
    }

    MouseArea {
        anchors.fill: parent
        onClicked: colorDialog.open()
    }

    ColorDialog {
        id: colorDialog
        selectedColor: root.color
        onAccepted: {
            root.color = selectedColor
            root.updated()
        }
    }
}
