import QtQuick
import QtQuick.Layouts

import SolTrace

STFormPanel {
    id: root

    property var module
    property bool singleColumn: false
    property int labelAlignment: Qt.AlignRight | Qt.AlignVCenter

    STPropertySeparator {
        title: "Placement"
    }

    STSwitch {
        text: "Use 3D widget"
        checked: App.view.mouse_mode === ViewModule.EditElement
        onToggled: App.view.mouse_mode = checked ? ViewModule.EditElement
                                                  : ViewModule.Camera
    }

    InlineDocumentation {
        key: "configure.layout.coordinates"
        Layout.fillWidth: true
    }

    STFormRow {
        label: "Coordinates"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STComboBox {
            id: positionModeCombo
            Layout.fillWidth: true
            model: ["Local", "Global"]
        }
    }

    StackLayout {
        id: positionSwipe
        Layout.fillWidth: true
        currentIndex: positionModeCombo.currentIndex

        ElementPositionEditor {
            Layout.preferredWidth: positionSwipe.width
            module: root.module
            singleColumn: root.singleColumn
            labelAlignment: root.labelAlignment
            targetProperty: "position"
            resetText: "Reset Local Default Position"
        }

        ElementPositionEditor {
            Layout.preferredWidth: positionSwipe.width
            module: root.module
            singleColumn: root.singleColumn
            labelAlignment: root.labelAlignment
            targetProperty: "global_position"
            resetText: "Reset Global Default Position"
        }
    }

    ElementRotationEditor {
        Layout.fillWidth: true
        module: root.module
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment
    }
}
