import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    Layout.fillHeight: true
    Layout.fillWidth: true
    spacing: 12

    Header {
        text: "Logs"
    }

    RowLayout {
        Layout.fillWidth: true

        STButton {
            text: "Open Log Directory"
            left_text_icon: "\uf07c"

            onClicked: AppData.log_list.open_log_directory()
        }

        Item {
            Layout.fillWidth: true
        }
    }

    ListView {
        id: log_list_view

        Layout.fillHeight: true
        Layout.fillWidth: true

        model: AppData.log_list
        clip: true
        spacing: 5

        ScrollBar.vertical: STScrollBar { }

        onCountChanged: log_list_view.positionViewAtEnd()

        delegate: Label {
            required property string content

            width: ListView.view.width
            wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            text: content
        }
    }
}
