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

            icon: "\uf00d"
            onClicked: {
                App.view.toggle_full_panel(root.available_width)
            }

            z: 3
        }

        StackLayout {
            currentIndex: App.view.full_panel.mode
            onCurrentIndexChanged: App.view.full_panel.mode = currentIndex

            anchors.fill: parent

            AdaptiveEditor {
                id: settingsView

                wideThreshold: 500
                listWidth: 200
                currentIndex: App.view.full_panel.settings_section
                onCurrentIndexChanged: {
                   App.view.full_panel.settings_section = currentIndex
                }

                model: ListModel {
                    ListElement { name: "Theme"; icon: "\uf53f" }
                    ListElement { name: "Scene"; icon: "\uf1b2" }
                    ListElement { name: "Language"; icon: "\uf0ac" }
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
                    currentIndex: settingsView.currentIndex

                    ThemeSettings {}

                    SceneSettings {}

                    LanguageSettings {}

                }
            }

            AdaptiveEditor {
                id: docsView

                wideThreshold: 500
                listWidth: 200
                currentIndex: App.view.full_panel.docs_section
                onCurrentIndexChanged: {
                    App.view.full_panel.docs_section = currentIndex
                }

                model: ListModel {
                    ListElement { name: "Overview"; icon: "\ue0bb" }
                    ListElement { name: "GitHub"; icon: "\uf02d" }
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
                    currentIndex: docsView.currentIndex

                    OverviewModule {}

                    GitHubModule {}

                    TeamSettings {}

                }
            }

            AdaptiveEditor {
                id: buildView

                wideThreshold: 500
                listWidth: 200
                currentIndex: App.view.full_panel.build_section
                onCurrentIndexChanged: {
                   App.view.full_panel.build_section = currentIndex
                }

                model: ListModel {
                    ListElement { name: "What's New"; icon: "✦" }
                    ListElement { name: "Build"; icon: "\uf7d9" }
                    ListElement { name: "Logs"; icon: "\uf15c" }
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
                    currentIndex: buildView.currentIndex

                    FeaturePanel {}

                    BInfoModule {}

                    LogSettings {}

                }
            }
        }


    }
}
