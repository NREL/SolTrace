import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

RowLayout {
    id: root

    spacing: 3

    signal itemSelected(db_entity entity);
    signal itemDeselected()

    STClickableLabel {
        text: qsTr("All")
        font.bold: true
        onClicked: {
            AppData.layout.clear_viewed_element()
            AppData.layout.clear_edited_element()
        }
    }

    // STIconButton {
    //     text: "\uf122"

    // }

    Repeater {
        model: AppData.layout.breadcrumb_model

        delegate: RowLayout {
            id: del_root
            required property string name
            required property var entity

            spacing: 3

            Label {
                id: label
                text: "\uf105"
                font.family: "Font Awesome 7 Free"

                font.pointSize: 16
            }

            STClickableLabel {
                text: del_root.name

                font.bold: true

                onClicked: {
                    AppData.layout.viewed_element = del_root.entity
                    AppData.layout.clear_edited_element()
                    root.itemSelected(del_root.entity)
                }
            }
        }
    }
}
