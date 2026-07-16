import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var current_db: AppData.file_source.current_database

    spacing: 8

    STPropertyPanel {
        Layout.fillWidth: true
        title: "Current Scene"

        RowLayout {
            Layout.fillWidth: true
            enabled: !!root.current_db

            STTextField {
                text: root.current_db ? root.current_db.name : ""
                Layout.fillWidth: true

                onTextEdited: {
                    if (root.current_db) {
                        root.current_db.name = text
                    }
                }
            }

            Label {
                font.family: "Font Awesome 7 Free"
                font.pointSize: 16
                text: "\uf303"
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

    STPropertyPanel {
        Layout.fillWidth: true
        Layout.fillHeight: true
        title: "Available Scenes"

        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true

            ListView {
                Layout.fillHeight: true
                Layout.fillWidth: true

                model: AppData.file_source
                clip: true

                ScrollBar.vertical: STScrollBar { }

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

            SceneListFooter {
                Layout.fillWidth: true
            }
        }
    }
}
