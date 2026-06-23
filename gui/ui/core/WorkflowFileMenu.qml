import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

STMenu {
    id: root

    FileController {
        id: file_controller
    }

    MenuItem {
        text: "New Scene"
        onTriggered: file_controller.load_new()
    }

    MenuItem {
        text: "Open File..."
        enabled: !AppData.file_source.is_loading
        onTriggered: file_controller.open_file()
    }

    STMenu {
        id: recents_menu

        title: "Recent Files"
        enabled: recent_instantiator.count > 0

        Instantiator {
            id: recent_instantiator
            model: file_controller.recent_files

            delegate: MenuItem {
                implicitWidth: 220
                text: file_controller.file_name(modelData)
                onTriggered: file_controller.load_recent(modelData)
            }

            onObjectAdded: (index, object) => recents_menu.insertItem(index, object)
            onObjectRemoved: (index, object) => recents_menu.removeItem(object)
        }
    }
}
