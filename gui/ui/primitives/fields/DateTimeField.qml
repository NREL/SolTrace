import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ColumnLayout {
    id: root
    property int month: 12
    property int day: 10
    property int hour: 12
    property int minute: 0
    signal modified()

    STPropertyPanel {
        id: date
        Layout.fillWidth: true

        title: qsTr("Date")
        collapsible: false

        GridLayout {
            rows: 2
            columns: 2
            Layout.fillWidth: true

            STSpinBoxField {
                id: monthField
                Layout.row: 0
                Layout.column: 0
                Layout.fillWidth: true
                Layout.maximumWidth: 200
                label: qsTr("Month")
                from: 1
                to: 12
                stepSize: 1
                onValueModified: { root.month = value; root.modified() }
            }
            Binding { monthField.value: root.month }

            STSpinBoxField {
                id: dayField
                Layout.row: App.view.left_panel.size == SplitPanelData.Small ? 1 : 0
                Layout.column: App.view.left_panel.size == SplitPanelData.Small ? 0 : 1
                Layout.preferredWidth: 200
                label: qsTr("Day")
                from: 1
                to: 31
                stepSize: 1
                onValueModified: { root.day = value; root.modified() }
            }
            Binding { dayField.value: root.day }
        }
    }

    STPropertyPanel {
        id: time
        Layout.fillWidth: true

        title: qsTr("Time")
        collapsible: false

        GridLayout {
            rows: 2
            columns: 2
            Layout.fillWidth: true

            STSpinBoxField {
                id: hourField
                Layout.row: 0
                Layout.column: 0
                Layout.fillWidth: true
                Layout.maximumWidth: 200
                label: qsTr("Hour")
                from: 0
                to: 23
                stepSize: 1
                onValueModified: { root.hour = value; root.modified() }
            }
            Binding { hourField.value: root.hour }

            STSpinBoxField {
                id: minuteField
                Layout.row: App.view.left_panel.size == SplitPanelData.Small ? 1 : 0
                Layout.column: App.view.left_panel.size == SplitPanelData.Small ? 0 : 1
                Layout.fillWidth: true
                Layout.maximumWidth: 200
                label: qsTr("Minute")
                from: 0
                to: 59
                stepSize: 1
                onValueModified: { root.minute = value; root.modified() }
            }
            Binding { minuteField.value: root.minute }
        }
    }
}
