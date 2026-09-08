import QtQuick.Layouts

import SolTrace

STFormPanel {
    id: root

    property var module
    property bool singleColumn: false

    ColorPickerField {
        id: elementColorPicker
        Layout.preferredWidth: 200
        color: root.module.color
        label: "Element Tint"
        onUpdated: {
            root.module.color = elementColorPicker.color
        }
    }

    STSwitch {
        text: "Show in 3D View"
        checked: !root.module.hidden
        onToggled: root.module.hidden = !checked
    }

    STSwitch {
        text: "Exclude from Simulation"
        checked: root.module.disabled
        onToggled: root.module.disabled = checked
    }

    STSwitch {
        text: "Ray Probe"
        checked: root.module.virtual_element
        onToggled: root.module.virtual_element = checked
    }
}
