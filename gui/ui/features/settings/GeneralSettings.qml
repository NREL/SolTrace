import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root

    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth
    clip: true

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Settings"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Documentation"
            collapsed: false

            ColumnLayout {
                Layout.fillWidth: true
                Layout.columnSpan: 2
                spacing: 8

                Label {
                    Layout.fillWidth: true
                    text: "Inline documentation shows short explanations next to controls throughout the app. It can be useful when learning SolTrace or when you want quick reminders about less familiar options."
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                }

                STSwitch {
                    text: "Show Inline Documentation"
                    checked: App.view.inline_docs
                    onToggled: App.view.inline_docs = checked
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
