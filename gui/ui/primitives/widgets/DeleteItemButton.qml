import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STIconButton {
    id: root
    icon: "\uf2ed"
    label: "Delete"

    onClicked: delete_item_pop.open()

    required property string title
    required property string itemType
    required property var replacementModel
    required property db_entity toDelete

    property bool isDangerous: false

    signal beforeOpened()
    signal deleteRequested()
    signal deleteReassignRequested(db_entity entity)

    STDialog {
        id: delete_item_pop

        title: parent.title

        onAccepted: {
            root.deleteRequested()
            delete_item_pop.close()
        }

        onAboutToShow: root.beforeOpened()

        ColumnLayout {
            anchors.fill: parent

            Label {
                visible: root.isDangerous
                Layout.fillWidth: true
                text: `This %1 is assigned to active users, how should these be handled?
                Delete and Unset will remove their %1 membership.
                Alternatively, you can delete this %1 and reassign those users to another %1.`
                .arg(root.itemType)
            }

            Label {
                visible: !root.isDangerous
                Layout.fillWidth: true
                text: "Confirm?"
            }
        }

        footer: STDialogButtonBox {
            STButton {
                text: qsTr("Cancel")
                DialogButtonBox.buttonRole: DialogButtonBox.RejectRole
            }

            STDangerousButton {
                visible: root.isDangerous
                text: qsTr("Reassign & Delete")
                DialogButtonBox.buttonRole: DialogButtonBox.DestructiveRole

                onClicked: select_replacement.open()

                SelectItemPopup {
                    id: select_replacement
                    source_model: root.replacementModel

                    exclude: [root.toDelete]

                    onSelectedEntity: function(entity) {
                        root.deleteReassignRequested(entity)
                        delete_item_pop.close()
                    }
                }
            }

            STButton {
                text: root.isDangerous ? qsTr("Unset & Delete") : qsTr("Delete")
                Material.foreground: Material.Red
                DialogButtonBox.buttonRole: DialogButtonBox.DestructiveRole
                idle_color: App.theme.destructiveGlassColor
                down_color: Material.color(Material.Red)

                onClicked: delete_item_pop.accept()
            }
        }
    }
}
