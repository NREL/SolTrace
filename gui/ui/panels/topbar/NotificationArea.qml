import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

STPopup {
    id: root

    property int new_notification_count: 0
    property int new_alert_count: 0

    onOpened: {
        new_notification_count = 0
        new_alert_count = 0
    }

    signal next_notification(message: string, type: int)

    ListModel {
        id: notification_model
    }

    Connections {
        target: AppData

        function onNotification(new_note) {

            if (!root.visible) {
                root.new_notification_count = root.new_notification_count + 1
                if (new_note.type > 0) {
                    root.new_alert_count = root.new_alert_count + 1
                }
            }

            console.log("Adding note: ", new_note.message)

            notification_model.insert(
                0,
                {
                    "message": new_note.message,
                    "type" : new_note.type,
                    "date" : Qt.formatDateTime(new Date(), "MMM d h:mm:ss AP"),
                }
            )

            root.next_notification(new_note.message, new_note.type)

            while (notification_model.count > 50) {
                notification_model.remove(notification_model.count - 1)
            }
        }
    }

    ColumnLayout {
        anchors.fill: parent

        ListView {
            Layout.preferredHeight: 300
            Layout.preferredWidth: 250
            Layout.fillWidth: true

            clip: true

            model: notification_model

            ScrollBar.vertical: STScrollBar { }

            spacing: 8

            delegate: GridLayout {
                id: note_delegate
                required property string message
                required property int type
                required property string date
                width: ListView.view.width

                property color message_color: {
                    switch (type) {
                    case 0: return Material.foreground
                    case 1: return Material.color(Material.Yellow)
                    case 2: return Material.color(Material.Red)
                    default: return Material.foreground
                    }
                }

                rowSpacing: 8
                columnSpacing: 8
                columns: 2

                Label {
                    fontSizeMode: Label.Fit
                    font.pointSize: 18
                    Layout.preferredWidth: 24
                    Layout.rowSpan: 2

                    horizontalAlignment: Qt.AlignHCenter
                    verticalAlignment: Qt.AlignVCenter

                    font.family: "Font Awesome 7 Free"

                    Material.foreground: note_delegate.message_color

                    text: {
                        switch (type) {
                        case 0: return "\uf05a"
                        case 1: return "\uf071"
                        case 2: return "\uf057"
                        default: return "\uf05a"
                        }
                    }
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                    elide: Label.ElideRight

                    Material.foreground: note_delegate.message_color

                    text: message
                }

                Label {
                    Layout.fillWidth: true
                    opacity: .75
                    text: note_delegate.date

                    Rectangle {
                        height: 1
                        color: Material.dividerColor
                        anchors.left: parent.left
                        anchors.right: parent.right
                        anchors.verticalCenter: parent.bottom
                        visible: notification_model.count > 1
                        anchors.verticalCenterOffset: 5
                    }
                }
            }

            Label {
                anchors.centerIn: parent

                enabled: false

                text: "No notifications"

                visible: parent.count === 0
            }
        }
    }

}
