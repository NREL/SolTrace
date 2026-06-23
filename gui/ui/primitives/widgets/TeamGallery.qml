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
            SubHeader { text: name }
            Label { text: role }
        }

        preview: ColumnLayout {
            property string name
            property string role
            property string description
            property string website
            property string email

            spacing: 8
            SubHeader { text: name }
            Label { text: role }
            Label {
                Layout.fillWidth: true
                text: description
                wrapMode: Text.WordWrap
            }
            RowLayout {
                STButton {
                    text: "Website"
                    text_icon: "\uf0c1"
                    onClicked: Qt.openUrlExternally(website)
                }
            }
        }
    }
}
