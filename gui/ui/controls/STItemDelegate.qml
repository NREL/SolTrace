import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ItemDelegate {
    id: control
    width: ListView.view ? ListView.view.width : implicitWidth
    highlighted: ListView.isCurrent ? ListView.isCurrent : false

    hoverEnabled: true

    background: Rectangle {
        implicitHeight: 24
        implicitWidth: 100
        opacity: control.enabled ? 1 : 0.3
        color: control.down
               ?
                   Material.rippleColor
                 :
                   control.highlighted
                   ?
                       Material.highlightedRippleColor
                     :
                       control.hovered
                       ?
                           Qt.alpha(Material.highlightedRippleColor, 0.35)
                         :
                           "transparent"
        radius: 14
    }
}
