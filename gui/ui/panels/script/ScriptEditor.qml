import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Popup {
    id: root
    parent: Overlay.overlay

    property var module: AppData.script

    height: parent.height * .9
    width: parent.width * .9

    anchors.centerIn: Overlay.overlay

    modal: true

    onOpened: text_edit_area.text = root.module.code

    ColumnLayout {
        anchors.fill: parent
        ScrollView {
            Layout.fillWidth: true
            Layout.fillHeight: true

            TextArea {
                id: text_edit_area


                onTextEdited: {
                    root.module.code = text
                }
            }
        }

        Label {
            id: status_label

            property bool has_parse_error: root.module.parse_errors.length

            text: "Parse error"

            visible: has_parse_error
        }

        Label {
            Layout.fillWidth: true

            visible: text.length

            text: root.module.parse_errors.join("\n")
        }
    }
}
