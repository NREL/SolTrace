import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Material
import SolTrace

STPropertyPanel {
    id: root
    property var model: null

    CardGallery {
        model: root.model
        Layout.fillWidth: true
        Layout.columnSpan: 2

        delegate: ColumnLayout {
            property string name
            property string role
            spacing: 5
            SubHeader {
                Layout.fillWidth: true
                text: name
                wrapMode: Text.WordWrap
            }
            Label {
                Layout.fillWidth: true
                text: role
                wrapMode: Text.WordWrap
            }
        }

        preview: ColumnLayout {
            property string name
            property string role
            property string description
            property string website
            property string email

            spacing: 8
            SubHeader {
                Layout.fillWidth: true
                text: name
                wrapMode: Text.WordWrap
            }
            Label {
                Layout.fillWidth: true
                text: role
                wrapMode: Text.WordWrap
            }
            Label {
                Layout.fillWidth: true
                text: description
                wrapMode: Text.WordWrap
            }
            RowLayout {
                STButton {
                    text: "Website"
                    left_text_icon: "\uf0c1"
                    onClicked: Qt.openUrlExternally(website)
                }
            }
        }
    }
}
