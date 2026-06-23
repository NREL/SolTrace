import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    function member(key) {
        return {
            name: App.docs.get("team." + key, "name"),
            role: App.docs.get("team." + key, "role"),
            description: App.docs.get("team." + key),
            website: App.docs.get("team." + key, "website"),
            email: App.docs.get("team." + key, "email")
        }
    }

    function addMembers(model, keys) {
        for (let key of keys) model.append(member(key))
    }

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "SolTrace Team"
        }

        TeamGallery {
            title: "Principal Investigator"
            Layout.fillWidth: true
            model: ListModel {
                Component.onCompleted: root.addMembers(this, ["w_hamilton"])
            }
        }

        TeamGallery {
            title: "Backend Team"
            Layout.fillWidth: true
            model: ListModel {
                Component.onCompleted: root.addMembers(this,
                    ["m_wagner", "t_brown", "j_maack", "l_fang", "n_edwards"])
            }
        }

        TeamGallery {
            title: "Frontend Team"
            Layout.fillWidth: true
            model: ListModel {
                Component.onCompleted: root.addMembers(this,
                    ["n_brunhart_lupo", "r_shantivong", "k_gruchalla"])
            }
        }

        Item { Layout.fillHeight: true }
    }
}
