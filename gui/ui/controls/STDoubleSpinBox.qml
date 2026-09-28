import QtQuick
import QtQuick.Controls.Material
import SolTrace

DoubleSpinBox {
    id: control
    editable: true

    from: -Infinity
    to: Infinity

    up.indicator: Rectangle {
        x: control.mirrored ? 0 : parent.width - width
        height: parent.height
        implicitWidth: 32
        implicitHeight: 32
        color: control.up.pressed ? Qt.rgba(1, 1, 1, 0.1) : "transparent"
        radius: height / 2
        opacity: (control.hovered || control.activeFocus) ? 1 : 0

        Label {
            text: "\u002b"
            font.family: "Font Awesome 7 Free"
            font.pointSize: control.font.pointSize
            anchors.centerIn: parent
        }
    }

    down.indicator: Rectangle {
        x: control.mirrored ? parent.width - width : 0
        height: parent.height
        implicitWidth: 32
        implicitHeight: 32
        color: control.down.pressed ? Qt.rgba(1, 1, 1, 0.1) : "transparent"
        radius: height / 2
        opacity: (control.hovered || control.activeFocus) ? 1 : 0

        Label {
            text: "\uf068"
            font.family: "Font Awesome 7 Free"
            font.pointSize: control.font.pointSize
            anchors.centerIn: parent
        }
    }

    background: WellRectangle {
        implicitWidth: 140
        implicitHeight: 32
        radius: height / 2

        function compute_border_color() {
            if (control.activeFocus) return control.Material.accentColor
            if (!enabled) return control.Material.hintTextColor
            return "black"
        }

        border.color: compute_border_color()
    }
}


