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
        text: "Build & Logs"
    }

    STPropertyPanel {
        Layout.fillWidth: true

        title: "Version"
        collapsed: false

        ColumnLayout {
            spacing: 6

            Label {
                Layout.fillWidth: true
                text: AppData.current_build_info
                font.bold: true
                wrapMode: Text.WordWrap
            }

            STClickableLabel {
                Layout.fillWidth: true
                text: "Copy build info to clipboard"
                color: Material.accentColor
                borderWidth: 0

                onClicked: AppData.copy_build_info_to_clipboard()
            }
        }
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
        id: logListView

        Layout.fillHeight: true
        Layout.fillWidth: true

        model: AppData.log_list
        clip: true
        spacing: 5

        ScrollBar.vertical: STScrollBar { }

        onCountChanged: logListView.positionViewAtEnd()

        delegate: Label {
            required property string content

            width: ListView.view.width
            wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            text: content
        }
    }
}
