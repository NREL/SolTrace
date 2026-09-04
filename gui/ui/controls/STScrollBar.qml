import QtQuick
import QtQuick.Controls

ScrollBar {
    id: control

    // ScrollIndicator is display-only. Desktop scroll surfaces use this
    // control so the thumb can be clicked and dragged and the track clicked.
    policy: ScrollBar.AsNeeded
    interactive: true
    hoverEnabled: true
    minimumSize: 0.08

    // Keep a comfortable pointer target even when the styled thumb is narrow.
    implicitWidth: orientation === Qt.Vertical ? 14 : 100
    implicitHeight: orientation === Qt.Horizontal ? 14 : 100
}
