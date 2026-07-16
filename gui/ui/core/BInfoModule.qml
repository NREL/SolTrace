import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Build Information"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Version"
            collapsed: false

            ColumnLayout {
                spacing: 6

                Label {
                    Layout.fillWidth: true
                    text: AppData.current_version_info
                    font.bold: true
                    elide: Text.ElideRight
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
    }
}

