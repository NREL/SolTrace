import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

Menu {
    background: Rectangle {
        implicitWidth: 150
        implicitHeight: 40

        radius: 14

        color: Qt.alpha(Material.backgroundColor, .90)
    }
}