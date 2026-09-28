import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

STPropertyPanel {
    id: root
    property var model: null

    Flow {
        Layout.fillWidth: true
        Layout.columnSpan: 2
        spacing: 8

        Repeater {
            model: root.model

            Rectangle {
                id: card

                required property string name
                required property string role
                required property string description
                required property string website
                required property string email

                width: 300
                height: cardContent.implicitHeight + 24
                color: App.theme.glassColorA(0.1)
                border.color: Qt.alpha("black", .25)
                radius: 8

                ColumnLayout {
                    id: cardContent
                    width: parent.width - 32
                    anchors.centerIn: parent
                    spacing: 5

                    SubHeader {
                        Layout.fillWidth: true
                        text: card.name
                        wrapMode: Text.WordWrap
                    }

                    Label {
                        Layout.fillWidth: true
                        text: card.role
                        wrapMode: Text.WordWrap
                    }
                }

                STIconButton {
                    anchors.top: parent.top
                    anchors.right: parent.right
                    anchors.margins: 8
                    icon: "\uf08e"
                    onClicked: detailsPopup.open()
                }

                Popup {
                    id: detailsPopup
                    width: 400
                    modal: true
                    anchors.centerIn: Overlay.overlay
                    padding: 12
                    closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

                    Overlay.modal: Rectangle {
                        color: Qt.rgba(0, 0, 0, 0.25)
                    }

                    ColumnLayout {
                        anchors.left: parent.left
                        anchors.right: parent.right
                        spacing: 8

                        SubHeader {
                            Layout.fillWidth: true
                            text: card.name
                            wrapMode: Text.WordWrap
                        }

                        Label {
                            Layout.fillWidth: true
                            text: card.role
                            wrapMode: Text.WordWrap
                        }

                        Label {
                            Layout.fillWidth: true
                            text: card.description
                            wrapMode: Text.WordWrap
                        }

                        RowLayout {
                            STButton {
                                text: "Website"
                                left_text_icon: "\uf0c1"
                                enabled: card.website.length > 0
                                onClicked: Qt.openUrlExternally(card.website)
                            }
                        }
                    }
                }
            }
        }
    }
}
