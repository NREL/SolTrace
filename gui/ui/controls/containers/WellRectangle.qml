import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Rectangle {
    id: root
    radius: height / 2
    color: App.theme.glassColor

    border.color: "black"
    border.width: 0
    
    Rectangle {
        anchors.fill: parent
        anchors.margins: 0
        
        radius: height / 2
        color: "transparent"
        border.color: Qt.alpha(root.border.color, .25)
        
        Rectangle {
            anchors.fill: parent
            anchors.margins: 0
            
            radius: height / 2
            color: "transparent"
            border.color: Qt.alpha(root.border.color, .1)
        }
    }
}
