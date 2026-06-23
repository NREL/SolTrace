import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

DialogButtonBox {
    visible: count > 0

    delegate: STButton { }

    //Material.foreground: parent.Material.foreground
    Material.foreground: undefined

    background: Item {
        implicitHeight: Material.dialogButtonBoxHeight

        Rectangle {
            color: Material.dividerColor
            anchors.left: parent.left
            anchors.leftMargin: 3
            anchors.right: parent.right
            anchors.rightMargin: 3
            height: 1
            anchors.top: parent.top
        }
    }
}