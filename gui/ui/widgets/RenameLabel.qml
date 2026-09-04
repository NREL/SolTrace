import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STClickableLabel {
    id: root
    onClicked: rename_pop.open()

    signal accepted(string name)

    STDialog {
        id: rename_pop

        modal: false

        STTextField {
            id: text_input
            placeholderText: qsTr("New name...")

            onAccepted: rename_pop.accept()
        }

        standardButtons: Dialog.Ok | Dialog.Cancel

        onAccepted: root.accepted(text_input.text)
    }

    Label {
        text: "\uf303"
        font.family: "Font Awesome 7 Free"

        font.pointSize: 16

        anchors.left: parent.right
        anchors.bottom: parent.top

        opacity: parent.containsMouse * .80

        Behavior on opacity {
            NumberAnimation {
                duration: 100
            }
        }
    }
}
