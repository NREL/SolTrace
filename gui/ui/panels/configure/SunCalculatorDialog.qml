import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STPopup {
    id: root

    modal: true

    width: 500

    property int calculator: 0
    property real latitude: 0
    property real longitude: 0
    property int year: 2026
    property int month: 1
    property int day: 1
    property int hour: 12
    property int minute: 0
    property int second: 0
    property int timezoneOffset: -7
    property real altitude: 1000
    property real pressure: 1013.25
    property real temperature: 20

    property string error_string: ""

    function openForCurrent() {
        calculator = App.sun.calc_data.calculator
        latitude = App.sun.calc_data.latitude
        longitude = App.sun.calc_data.longitude
        year = App.sun.calc_data.year
        month = App.sun.calc_data.month
        day = App.sun.calc_data.day
        hour = App.sun.calc_data.hour
        minute = App.sun.calc_data.minute
        second = App.sun.calc_data.second
        timezoneOffset = App.sun.calc_data.timezone_offset
        altitude = App.sun.calc_data.altitude
        pressure = App.sun.calc_data.pressure
        temperature = App.sun.calc_data.temperature
        open()
    }

    function apply_changes() {
        root.error_string = App.sun.apply_calculator(calculator,
                                              latitude,
                                              longitude,
                                              year,
                                              month,
                                              day,
                                              hour,
                                              minute,
                                              second,
                                              timezoneOffset,
                                              altitude,
                                              pressure,
                                              temperature)

        if (root.error_string.length > 0) {
            return;
        }

        root.close()
    }

    contentItem: ColumnLayout {
        spacing: 10
        //implicitWidth: Math.min(root.parent ? root.parent.width - 36 : 420, 420)

        RowLayout {
            Layout.fillWidth: true

            Label {
                text: "Calculator"
            }

            STComboBox {
                Layout.fillWidth: true
                currentIndex: root.calculator
                model: ["Legacy", "Duffie", "SOLPOS", "SPA"]
                onCurrentIndexChanged: {
                    root.calculator = currentIndex
                }
            }
        }

        GridLayout {
            Layout.fillWidth: true
            columns: 2

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Latitude"
                from: -90
                to: 90
                decimals: 4
                suffix: "deg"
                value: root.latitude
                onValueModified: root.latitude = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Longitude"
                from: -180
                to: 180
                decimals: 4
                suffix: "deg"
                value: root.longitude
                onValueModified: root.longitude = value
            }

            Rectangle {
                Layout.columnSpan: 2
                Layout.fillWidth: true
                height: 1
                color: Material.dividerColor
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Year"
                from: 1900
                to: 2200
                value: root.year
                onValueModified: root.year = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Month"
                from: 1
                to: 12
                value: root.month
                onValueModified: root.month = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Day"
                from: 1
                to: 31
                value: root.day
                onValueModified: root.day = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Hour"
                from: 0
                to: 23
                value: root.hour
                onValueModified: root.hour = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Minute"
                from: 0
                to: 59
                value: root.minute
                onValueModified: root.minute = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Second"
                from: 0
                to: 59
                value: root.second
                onValueModified: root.second = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "UTC Offset"
                from: -12
                to: 14
                value: root.timezoneOffset
                onValueModified: root.timezoneOffset = value
            }

            Item {}

            Rectangle {
                Layout.columnSpan: 2
                Layout.fillWidth: true
                height: 1
                color: Material.dividerColor
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Altitude"
                from: -440
                to: 8850
                decimals: 1
                suffix: "m"
                value: root.altitude
                onValueModified: root.altitude = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Pressure"
                from: 950
                to: 1050
                decimals: 2
                value: root.pressure
                onValueModified: root.pressure = value
            }

            STSpinBoxField {
                Layout.fillWidth: true
                label: "Temperature"
                from: -50
                to: 60
                decimals: 1
                value: root.temperature
                onValueModified: root.temperature = value
            }
        }

        Label {
            text: root.error_string
            visible: text.length > 0
            Layout.fillWidth: true
            wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            color: Material.color(Material.Yellow)
        }

        RowLayout {
            STButton {
                text: "Apply"
                onClicked: root.apply_changes()
            }
        }
    }
}
