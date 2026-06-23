import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var current_db: AppData.file_source.current_database

    spacing: 8

    Label {
        Layout.fillWidth: true

        text: "Available Scenes"
        font.bold: true
    }

    Rectangle {
        color: Material.dividerColor
        Layout.preferredHeight: 1
        Layout.fillWidth: true
        Layout.leftMargin: 3
        Layout.rightMargin: 3
    }

    ListView {
        Layout.fillHeight: true
        Layout.fillWidth: true

        model: AppData.file_source
        clip: true

        delegate: STItemDelegate {
            id: db_delegate

            required property int index
            required property var database

            text: database.name
            font.pointSize: 18
            highlighted: root.current_db === database

            onClicked: {
                AppData.file_source.set_current(index)
                App.view.simulation_content_view = false
            }
        }
    }

    Rectangle {
        color: Material.dividerColor
        Layout.preferredHeight: 1
        Layout.fillWidth: true
        Layout.leftMargin: 3
        Layout.rightMargin: 3
    }

    Label {
        text: "Current Scene"
        font.bold: true
    }

    RowLayout {
        enabled: !!root.current_db

        Label {
            font.family: "Font Awesome 7 Free"
            font.pointSize: 16
            text: "\uf303"
        }

        STTextField {
            text: root.current_db ? root.current_db.name : ""
            Layout.fillWidth: true

            onTextChanged: {
                if (root.current_db) {
                    root.current_db.name = text
                }
            }
        }
    }

    SceneListFooter {
        Layout.fillWidth: true
    }
}
