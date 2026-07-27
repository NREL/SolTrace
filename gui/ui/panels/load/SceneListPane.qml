import QtCore
import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Dialogs
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var current_db: AppData.file_source.current_database

    spacing: 8

    STPropertyPanel {
        Layout.fillWidth: true
        title: "Current Scene"

        enabled: !!root.current_db

        RowLayout {
            Layout.fillWidth: true
            Layout.columnSpan: 2

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

        RowLayout {
            Layout.fillWidth: true
            Layout.columnSpan: 2

            Item {
                Layout.fillWidth: true
            }


            STIconButton {
                enabled: !AppData.file_source.is_loading
                icon: "\uf56e"
                label: "Export"
                toolTip: "Export scene to file"

                onClicked: {
                    if (!AppData.file_source.save_current_dialog()) {
                        file_dialog.open()
                    }
                }

                FileDialog {
                    id: file_dialog

                    fileMode: FileDialog.SaveFile

                    defaultSuffix: "json"

                    nameFilters: ["JSON files (*.json)"]

                    onAccepted:
                        AppData.file_source.save_current(file_dialog.selectedFile)


                    Settings {
                        id: save_file_history

                        category: "save_file_history"

                        property alias last_selected_file: file_dialog.selectedFile
                        property alias last_selected_folder: file_dialog.currentFolder
                    }
                }
            }
        }
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
