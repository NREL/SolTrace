import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Controls.Material
import SolTrace

ShadowedGlassRectangle {
    id: root
    required property int available_width

    Item {

        anchors.fill: parent
        anchors.margins: 8

        STIconButton {
            id: close_button

            anchors.top: parent.top
            anchors.right: parent.right

            text: "\uf00d"
            onClicked: {
                App.view.settings_panel.visible = false
                App.view.left_panel.visible = App.view.left_panel.saved_visible
                App.view.right_panel.visible = App.view.right_panel.saved_visible
            }

            z: 3
        }

        AdaptiveEditor {
            id: editor

            anchors.fill: parent

            wideThreshold: 500
            listWidth: 200
            currentIndex: App.view.documentation_section
            onCurrentIndexChanged: {
                App.view.documentation_section = currentIndex
            }

            model: ListModel {
                ListElement { name: "What's New"; icon: "✦" }
                ListElement { name: "Theme"; icon: "\uf53f" }
                ListElement { name: "Scene"; icon: "\uf1b2" }
                ListElement { name: "Logs"; icon: "\uf15c" }
                ListElement { name: "Documentation"; icon: "\uf02d" }
                ListElement { name: "Team"; icon: "\uf500" }
            }

            listDelegate: ItemDelegate {
                text: itemModel ? itemModel.name : ""
                highlighted: isCurrent
                width: parent ? parent.width : implicitWidth

                contentItem: RowLayout {
                    spacing: 8
                    Label {
                        text: itemModel ? itemModel.icon : ""
                        font.family: "Font Awesome 7 Free"
                        font.pointSize: 14
                    }
                    Label {
                        text: itemModel ? itemModel.name : ""
                        Layout.fillWidth: true
                    }
                }

                background: Rectangle {
                    implicitHeight: 36
                    implicitWidth: 100
                    opacity: enabled ? 1 : 0.3
                    color: parent.down ? Material.rippleColor
                         : parent.highlighted ? Qt.rgba(Material.accentColor.r,
                                                         Material.accentColor.g,
                                                         Material.accentColor.b, 0.12)
                         : "transparent"
                    radius: 14
                }
            }

            detailView: StackLayout {
                currentIndex: editor.currentIndex

                FeaturePanel {}

                ThemeSettings {}

                SceneSettings {}

                LogSettings {}

                DocumentationSettings {}

                TeamSettings {}

            }
        }
    }
}
