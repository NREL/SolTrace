import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

STPropertyPanel {
    id: root
    property real latitude: 0
    property real longitude: 0
    signal modified()

    Layout.fillWidth: true

    title: "Location"
    collapsible: false

    GridLayout {
        rows: 2
        columns: 2
        Layout.fillWidth: true

        STSpinBoxField {
            id: latField
            Layout.row: 0
            Layout.column: 0
            Layout.fillWidth: true
            Layout.maximumWidth: 200
            label: "Latitude"
            from: -90
            to: 90
            decimals: 4
            suffix: "deg"
            onValueModified: { root.latitude = value; root.modified() }
        }
        Binding { latField.value: root.latitude }

        STSpinBoxField {
            id: lonField
            Layout.row: App.view.left_panel.size == PanelData.Small ? 1 : 0
            Layout.column: App.view.left_panel.size == PanelData.Small ? 0 : 1
            Layout.fillWidth: true
            Layout.maximumWidth: 200
            label: "Longitude"
            from: -180
            to: 180
            decimals: 4
            suffix: "deg"
            onValueModified: { root.longitude = value; root.modified() }
        }
        Binding { lonField.value: root.longitude }
    }
}
