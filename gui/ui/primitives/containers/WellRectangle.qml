import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Rectangle {
    radius: height / 2
    color: App.theme.glassColor
    
    Rectangle {
        anchors.fill: parent
        anchors.margins: 0
        
        radius: height / 2
        color: "transparent"
        border.color: Qt.alpha("black", .25)
        
        Rectangle {
            anchors.fill: parent
            anchors.margins: 0
            
            radius: height / 2
            color: "transparent"
            border.color: Qt.alpha("black", .1)
        }
    }
}
