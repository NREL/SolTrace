import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts
import QtQuick.Effects

Dialog {
    id: root

    header: Label {
        text: root.title

        visible: parent?.parent === Overlay.overlay && root.title

        horizontalAlignment: Qt.AlignHCenter

        padding: 8
        bottomPadding: 4
        font.pixelSize: 16

        Rectangle {
            color: Material.dividerColor
            anchors.left: parent.left
            anchors.leftMargin: 3
            anchors.right: parent.right
            anchors.rightMargin: 3
            height: 1
            anchors.bottom: parent.bottom
        }
    }

    background: Rectangle {
        radius: 14

        color: Qt.alpha(Material.backgroundColor, .90)

        RectangularShadow {
            anchors.fill: parent
            offset.x: 0
            offset.y: 0
            radius: parent.radius
            blur: 30
            spread: 2
            color: Material.dropShadowColor
        }
    }

    footer: STDialogButtonBox { }

    Overlay.modal: Rectangle {
        //color: "#66000000"
        color: Material.backgroundDimColor
    }

    margins: 1
    padding: 8
}
