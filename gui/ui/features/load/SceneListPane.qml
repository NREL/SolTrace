import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var current_db: AppData.file_source.current_database
    property bool showFooter: true

    spacing: 8

    FileController {
        id: file_controller
    }

    Label {
        Layout.fillWidth: true
        text: "Current Scene"
        font.bold: true
    }

    RowLayout {
        Layout.fillWidth: true
        enabled: !!root.current_db
        spacing: 8

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            font.pointSize: 16
            text: "\uf303"
        }

        STTextField {
            Layout.fillWidth: true
            text: root.current_db ? root.current_db.name : ""

            onTextEdited: {
                if (root.current_db) {
                    root.current_db.name = text
                }
            }
        }

        STIconButton {
            enabled: !!root.current_db && !AppData.file_source.is_loading
            icon: "\uf0c7"
            label: "Save"
            toolTip: "Save scene to file"

            onClicked: file_controller.save_current()
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
        Layout.fillWidth: true
        text: "Available Scenes"
        font.bold: true
    }

    Label {
        Layout.fillWidth: true
        visible: scene_list.count === 0
        text: "No scenes loaded."
        opacity: 0.7
        horizontalAlignment: Text.AlignHCenter
    }

    ListView {
        id: scene_list

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

    Rectangle {
        color: Material.dividerColor
        Layout.preferredHeight: 1
        Layout.fillWidth: true
        Layout.leftMargin: 3
        Layout.rightMargin: 3
    }

    SceneListFooter {
        Layout.fillWidth: true
        visible: root.showFooter
    }
}
