import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

RowLayout {
    id: root

    required property var blur_source
    required property int available_width

    Layout.preferredWidth: implicitWidth

    STIconButton {
        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitHeight
        label.font.pointSize: 20

        text: "\uf0f3"
        toolTip: "Notifications"

        onClicked: notification_settings.open()

        Rectangle {
            id: notification_pill
            width: 8
            height: width
            radius: width / 2

            anchors.horizontalCenter: parent.left
            anchors.verticalCenter: parent.top

            opacity: notification_settings.new_alert_count > 0

            Connections {
                target: notification_settings
                function onNew_alert_countChanged() {
                    pulse_animation.restart()
                }
            }

            NumberAnimation {
                id: pulse_animation

                target: notification_pill
                property: "width"
                from: 16
                to: 8

                duration: 1000

                easing: Easing.OutElastic
            }

            color: Material.color(Material.Red)
        }

        NotificationArea {
            id: notification_settings

            onNext_notification: function(message, type) {
                new_note_pop.type = type
                new_note_pop.message = message
                new_note_pop.show_popup()
            }
        }

        STPopup {
            id: new_note_pop

            property string message
            property int type

            property color message_color: {
                switch (type) {
                case 0: return Material.foreground
                case 1: return Material.color(Material.Yellow)
                case 2: return Material.color(Material.Red)
                default: return Material.foreground
                }
            }

            contentWidth: 250
            contentHeight: 50

            function show_popup() {
                open()
                close_note_timer.restart()
            }

            closePolicy: Popup.NoAutoClose

            Timer {
                id: close_note_timer

                interval: 3000

                onTriggered: new_note_pop.close()
            }

            RowLayout {
                id: note_layout
                anchors.fill: parent

                Label {
                    Material.foreground: new_note_pop.message_color

                    font.pointSize: 18
                    Layout.preferredWidth: 24

                    horizontalAlignment: Qt.AlignHCenter
                    verticalAlignment: Qt.AlignVCenter

                    font.family: "Font Awesome 7 Free"

                    text: {
                        switch (new_note_pop.type) {
                        case 0: return "\uf05a"
                        case 1: return "\uf071"
                        case 2: return "\uf057"
                        default: return "\uf05a"
                        }
                    }
                }

                Label {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                    elide: Label.ElideRight

                    Material.foreground: new_note_pop.message_color

                    text: new_note_pop.message
                }
            }

            MouseArea {
                anchors.fill: note_layout
                onClicked: {
                    new_note_pop.close()
                    notification_settings.open()
                }
            }
        }
    }

    STIconButton {
        id: settings_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitHeight
        label.font.pointSize: 20

        text: "\uf013"
        toolTip: "Settings"

        Connections {
            target: App.view.settings_panel
            function onVisible_changed() {
                if (App.view.settings_panel.visible) {
                    App.view.left_panel.saved_visible = App.view.left_panel.visible
                    App.view.left_panel.visible = false

                    App.view.right_panel.saved_visible = App.view.right_panel.visible
                    App.view.right_panel.visible = false
                } else {
                    App.view.left_panel.visible = App.view.left_panel.saved_visible
                    App.view.right_panel.visible = App.view.right_panel.saved_visible
                    App.view.fit_panels(root.available_width, false, true)
                }
            }
        }

        onClicked: {
            App.view.settings_panel.visible = !App.view.settings_panel.visible
        }
    }

    STIconButton {
        id: rightpanel_open
        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitHeight

        text: "\uf0c9"
        toolTip: (App.view.right_panel.visible ? "Close": "Open") + " Right Panel"
        label.font.pointSize: 20
        onClicked: {
            if (App.view.settings_panel.visible) {
                App.view.settings_panel.visible = false
                return
            }

            App.view.right_panel.visible = !App.view.right_panel.visible && !App.view.settings_panel.visible
            App.view.fit_panels(root.available_width, true)
        }
    }
}
