import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STIconButton {
    id: root
    icon: "\uf055"
    label: qsTr("New")

    onClicked: new_item_pop.open()

    required property string title

    signal createRequested(string name)

    STDialog {
        id: new_item_pop

        title: parent.title

        onAccepted: {
            root.createRequested(new_name_fld.text)
            new_item_pop.close()
        }

        GridLayout {
            anchors.fill: parent

            columns: 2

            Label {
                text: qsTr("Name:")
            }

            STTextField {
                Layout.fillWidth: true
                id: new_name_fld

                onAccepted: new_item_pop.accept()
            }
        }

        standardButtons: Dialog.Ok | Dialog.Cancel
    }
}
