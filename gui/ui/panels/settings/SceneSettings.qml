import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root

    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth
    clip: true

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Scene"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Geometry"
            collapsed: false

            ColorPickerField {
                id: geometryColorPicker
                color: App.view.sim.geometry_color
                label: "Global Geometry Color"
                onUpdated: {
                    App.view.sim.geometry_color = geometryColorPicker.color
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
