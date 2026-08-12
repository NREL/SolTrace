import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    spacing: 8

    FileController {
        id: file_controller
    }

    RowLayout {
        Layout.fillWidth: true
        spacing: 8

        STButton {
            Layout.fillWidth: true
            enabled: !AppData.file_source.is_loading
            text: "Import"
            left_text_icon: "\uf56f"
            onClicked: file_menu.open()

            WorkflowFileMenu {
                id: file_menu
                showNewScene: false
            }
        }

        STButton {
            Layout.fillWidth: true
            text: "Examples"
            left_text_icon: "\ue52f"
            onClicked: file_controller.open_example()
        }

        STButton {
            Layout.fillWidth: true
            text: "Blank"
            left_text_icon: "\uf055"
            onClicked: new_name_pop.open()

            STDialog {
                id: new_name_pop

                modal: false

                STTextField {
                    id: text_input
                    placeholderText: "New scene name..."

                    onAccepted: new_name_pop.accept()
                }

                standardButtons: Dialog.Ok | Dialog.Cancel

                onAccepted: AppData.file_source.append_new(text_input.text)
            }
        }
    }

    SceneListPane {
        Layout.fillWidth: true
        Layout.fillHeight: true
        showFooter: false
    }
}
