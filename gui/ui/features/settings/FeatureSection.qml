import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STPropertyPanel {
    id: root

    property var keys: []

    collapsed: false

    ColumnLayout {
        Layout.fillWidth: true
        Layout.columnSpan: 2
        spacing: 8

        Repeater {
            model: root.keys

            ColumnLayout {
                required property string modelData

                Layout.fillWidth: true
                spacing: 4

                readonly property string docKey: "settings.features." + modelData

                SubHeader {
                    Layout.fillWidth: true
                    text: App.docs.get(docKey, "title")
                    wrapMode: Text.WordWrap
                }

                Label {
                    Layout.fillWidth: true
                    text: App.docs.get(docKey)
                    wrapMode: Text.WordWrap
                }
            }
        }
    }
}
