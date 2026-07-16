import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    STPropertyPanel {
        Layout.fillWidth: true

        title: "Geometry"
        collapsible: false

        ColorPickerField {
            id: geometryColorPicker
            color: App.view.sim.geometry_color
            label: "Global Geometry Color"
            onUpdated: {
                App.view.sim.geometry_color = geometryColorPicker.color
            }
        }
    }

}
