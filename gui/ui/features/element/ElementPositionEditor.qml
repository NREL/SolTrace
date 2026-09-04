import QtQuick
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var module
    property string targetProperty: "position"
    property string resetText: "Reset"
    property bool singleColumn: false
    property int labelAlignment: Qt.AlignRight | Qt.AlignVCenter

    property real xPosition: 0
    property real yPosition: 0
    property real zPosition: 0

    function set_position_values(x, y, z) {
        xPosition = x
        yPosition = y
        zPosition = z
        x_pos.value = x
        y_pos.value = y
        z_pos.value = z
    }

    function current_position() {
        return targetProperty === "global_position"
                ? module.global_position
                : module.position
    }

    function refresh_position() {
        const position = current_position()
        set_position_values(position.x, position.y, position.z)
    }

    function update_position() {
        const position = Qt.vector3d(xPosition, yPosition, zPosition)
        if (targetProperty === "global_position") {
            module.global_position = position
        } else {
            module.position = position
        }
    }

    function reset_position() {
        const position = Qt.vector3d(0, 0, 0)
        if (targetProperty === "global_position") {
            module.global_position = position
        } else {
            module.position = position
        }
        refresh_position()
    }

    Component.onCompleted: refresh_position()

    Connections {
        target: root.module

        function onPosition_changed() {
            if (root.targetProperty === "position")
                root.refresh_position()
        }

        function onGlobal_position_changed() {
            if (root.targetProperty === "global_position")
                root.refresh_position()
        }
    }

    STFormRow {
        label: "X"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STDoubleSpinBox {
            id: x_pos
            Layout.fillWidth: true
            from: -Infinity
            to: Infinity
            onValueModified: {
                root.xPosition = value
                root.update_position()
            }
            decimals: 4
        }
    }

    STFormRow {
        label: "Y"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STDoubleSpinBox {
            id: y_pos
            Layout.fillWidth: true
            from: -Infinity
            to: Infinity
            onValueModified: {
                root.yPosition = value
                root.update_position()
            }
            decimals: 4
        }
    }

    STFormRow {
        label: "Z"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STDoubleSpinBox {
            id: z_pos
            Layout.fillWidth: true
            from: -Infinity
            to: Infinity
            onValueModified: {
                root.zPosition = value
                root.update_position()
            }
            decimals: 4
        }
    }

    STButton {
        Layout.fillWidth: true
        text: root.resetText
        onClicked: root.reset_position()
    }
}
