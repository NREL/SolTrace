import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property string error_string: ""

    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    readonly property int fieldSpacing: 8

    function apply_changes() {
        root.error_string = App.sun.apply_calculator(
            App.sun.calc_data.calculator,
            App.sun.calc_data.latitude,
            App.sun.calc_data.longitude,
            App.sun.calc_data.year,
            App.sun.calc_data.month,
            App.sun.calc_data.day,
            App.sun.calc_data.hour,
            App.sun.calc_data.minute,
            App.sun.calc_data.second,
            App.sun.calc_data.timezone_offset,
            App.sun.calc_data.altitude,
            App.sun.calc_data.pressure,
            App.sun.calc_data.temperature)

        if (root.error_string.length > 0) {
            return;
        }
    }

    spacing: 10

    RowLayout {
        Layout.fillWidth: true

        Label {
            text: "Calculator"
        }

        STComboBox {
            Layout.fillWidth: true
            currentIndex: App.sun.calc_data.calculator
            model: ["Legacy", "Duffie", "SOLPOS", "SPA"]
            onCurrentIndexChanged: {
                App.sun.calc_data.calculator = currentIndex
            }
        }
    }

    Rectangle {
        Layout.fillWidth: true
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    Flow {
        id: positionGroup
        spacing: root.fieldSpacing
        Layout.fillWidth: true

        readonly property int fieldWidth: root.singleColumn ? parent.width : (parent.width - root.fieldSpacing) / 2

        STSpinBoxField {
            width: positionGroup.fieldWidth
            label: "Latitude"
            from: -90
            to: 90
            decimals: 4
            suffix: "deg"
            value: App.sun.calc_data.latitude
            onValueModified: App.sun.calc_data.latitude = value
        }

        STSpinBoxField {
            width: positionGroup.fieldWidth
            label: "Longitude"
            from: -180
            to: 180
            decimals: 4
            suffix: "deg"
            value: App.sun.calc_data.longitude
            onValueModified: App.sun.calc_data.longitude = value
        }
    }

    Rectangle {
        Layout.fillWidth: true
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    Flow {
        id: dateGroup
        spacing: root.fieldSpacing
        Layout.fillWidth: true

        readonly property int fieldWidth: root.singleColumn ? parent.width : (parent.width - 2 * root.fieldSpacing) / 3

        STSpinBoxField {
            width: dateGroup.fieldWidth
            label: "Year"
            from: 1900
            to: 2200
            value: App.sun.calc_data.year
            onValueModified: App.sun.calc_data.year = value
        }

        STSpinBoxField {
            width: dateGroup.fieldWidth
            label: "Month"
            from: 1
            to: 12
            value: App.sun.calc_data.month
            onValueModified: App.sun.calc_data.month = value
        }

        STSpinBoxField {
            width: dateGroup.fieldWidth
            label: "Day"
            from: 1
            to: 31
            value: App.sun.calc_data.day
            onValueModified: App.sun.calc_data.day = value
        }
    }

    Flow {
        id: timeGroup
        spacing: root.fieldSpacing
        Layout.fillWidth: true

        readonly property int fieldWidth: root.singleColumn ? parent.width : (parent.width - 100 - 3 * root.fieldSpacing) / 3

        STSpinBoxField {
            width: timeGroup.fieldWidth
            label: "Hour"
            from: 0
            to: 23
            value: App.sun.calc_data.hour
            onValueModified: App.sun.calc_data.hour = value
        }

        STSpinBoxField {
            width: timeGroup.fieldWidth
            label: "Minute"
            from: 0
            to: 59
            value: App.sun.calc_data.minute
            onValueModified: App.sun.calc_data.minute = value
        }

        STSpinBoxField {
            width: timeGroup.fieldWidth
            label: "Second"
            from: 0
            to: 59
            value: App.sun.calc_data.second
            onValueModified: App.sun.calc_data.second = value
        }

        STSpinBoxField {
            width: root.singleColumn ? parent.width : 100
            label: "UTC Offset"
            from: -12
            to: 14
            value: App.sun.calc_data.timezone_offset
            onValueModified: App.sun.calc_data.timezone_offset = value
        }
    }

    RowLayout {
        Layout.fillWidth: true

        Label {
            text: "Preset Season"
        }

        STComboBox {
            Layout.fillWidth: true
            currentIndex: -1
            model: ["Sprint", "Summer", "Fall", "Winter"]
            onCurrentIndexChanged: {
                if (currentIndex >= 0) {
                    [
                        App.sun.calc_data.set_spring,
                        App.sun.calc_data.set_summer,
                        App.sun.calc_data.set_fall,
                        App.sun.calc_data.set_winter,
                    ][currentIndex]()
                }
            }
        }
    }

    RowLayout {
        Layout.fillWidth: true

        Label {
            text: "Preset Time"
        }

        STComboBox {
            Layout.fillWidth: true
            currentIndex: -1
            model: ["Morning", "Noon", "Afternoon"]
            onCurrentIndexChanged: {
                if (currentIndex >= 0) {
                    [
                        App.sun.calc_data.set_morning,
                        App.sun.calc_data.set_noon,
                        App.sun.calc_data.set_afternoon,
                    ][currentIndex]()
                }
            }
        }
    }

    Rectangle {
        Layout.columnSpan: 2
        Layout.fillWidth: true
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    Flow {
        id: extraOptions
        spacing: root.fieldSpacing
        Layout.fillWidth: true

        visible: App.sun.calc_data.calculator === SolarCalculatorData.SPA

        readonly property int fieldWidth: root.singleColumn ? parent.width : (parent.width - 2 * root.fieldSpacing) / 3

        STSpinBoxField {
            width: extraOptions.fieldWidth
            label: "Altitude"
            from: -440
            to: 8850
            decimals: 1
            suffix: "m"
            value: App.sun.calc_data.altitude
            onValueModified: App.sun.calc_data.altitude = value
        }

        STSpinBoxField {
            width: extraOptions.fieldWidth
            label: "Pressure"
            from: 950
            to: 1050
            decimals: 2
            value: App.sun.calc_data.pressure
            onValueModified: App.sun.calc_data.pressure = value
        }

        STSpinBoxField {
            width: extraOptions.fieldWidth
            label: "Temperature"
            from: -50
            to: 60
            decimals: 1
            value: App.sun.calc_data.temperature
            onValueModified: App.sun.calc_data.temperature = value
        }
    }

    Label {
        text: root.error_string
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Label.WrapAtWordBoundaryOrAnywhere
        color: Material.color(Material.Yellow)
    }

    STButton {
        Layout.fillWidth: true
        text: "Apply"
        onClicked: root.apply_changes()
    }
}
