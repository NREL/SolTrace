import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

import SolTrace

RowLayout {
    id: root

    FileController {
        id: file_controller
    }

    STIconButton {
        enabled: !!AppData.file_source.current_database
        icon: "\uf2ed"
        toolTip: "Delete Scene"
        onClicked: AppData.file_source.delete_current()
    }

    Item {
        Layout.fillWidth: true
    }

    STIconButton {
        enabled: !!AppData.file_source.current_database
        icon: "\ue52f"
        label: "Example"
        toolTip: "Load Example"
        onClicked: file_controller.open_example()
    }

    STIconButton {
        enabled: !AppData.file_source.is_loading
        icon: "\uf56f"
        label: "Import"
        toolTip: "Import Scene"

        onClicked: file_menu.open()

        WorkflowFileMenu {
            id: file_menu
        }
    }

    STIconButton {
        icon: "\uf055"
        label: "Blank"
        toolTip: "Create Blank Scene"

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
