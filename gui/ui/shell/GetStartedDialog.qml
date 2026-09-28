import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STDialog {
    id: root

    title: "Get Started"

    modal: false

    closePolicy: Dialog.CloseOnEscape | Dialog.CloseOnPressOutside

    standardButtons: Dialog.Close

    parent: Overlay.overlay

    width: Math.min(720, Math.max(320, Overlay.overlay.width * 0.72))
    height: Math.min(620, Math.max(360, Overlay.overlay.height * 0.78))

    anchors.centerIn: parent

    ColumnLayout {
        anchors.fill: parent
        spacing: 8

        GetStartedModule {
            Layout.fillWidth: true
            Layout.fillHeight: true
        }

        Rectangle {
            color: Material.dividerColor
            height: 1
            Layout.fillWidth: true
        }

        ColumnLayout {
            Label {
                Layout.fillWidth: true
                text: "SolTrace has an in-line documentation system to provide useful information next to controls. If you are new to CSP software or SolTrace, you may consider enabling it."
                wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            }

            STSwitch {
                text: "Show Inline Documentation"
                checked: App.view.inline_docs
                onToggled: App.view.inline_docs = checked
            }
        }


        STSwitch {
            Layout.fillWidth: true
            text: "Show Get Started at startup"
            checked: AppData.show_get_started_on_start
            onToggled: AppData.show_get_started_on_start = checked
        }
    }
}
