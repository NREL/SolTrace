import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Documentation"
        }

        Rectangle {
            id: banner
            Layout.fillWidth: true
            implicitHeight: bannerContent.implicitHeight + 24
            radius: 12
            clip: true
            gradient: Gradient {
                orientation: Gradient.Horizontal
                GradientStop { position: 0.0; color: "#1a3a8a" }
                GradientStop { position: 0.35; color: "#1e2a6e" }
                GradientStop { position: 0.7; color: "#2d1b69" }
                GradientStop { position: 1.0; color: "#3b0764" }
            }

            ColumnLayout {
                id: bannerContent
                anchors.fill: parent
                anchors.margins: 12
                spacing: 8

                Label {
                    Layout.alignment: Qt.AlignHCenter
                    text: "Join SolTrace Discussions →"
                    color: "white"
                    font.bold: true
                    font.pointSize: App.theme.labelSize + 2

                    MouseArea {
                        anchors.fill: parent
                        cursorShape: Qt.PointingHandCursor
                        hoverEnabled: true
                        onClicked: Qt.openUrlExternally("https://github.com/NatLabRockies/SolTrace/discussions")
                        onContainsMouseChanged: parent.opacity = containsMouse ? 0.85 : 1.0
                    }
                }

                RowLayout {
                    id: bannerPills

                    Layout.preferredWidth: implicitWidth
                    Layout.alignment: Qt.AlignHCenter
                    spacing: 6

                    Repeater {
                        model: [
                            { emoji: "\uf0a1", label: "Announcements", slug: "announcements" },
                            { emoji: "\uf4ad", label: "General", slug: "general" },
                            { emoji: "\uf0eb", label: "Ideas", slug: "ideas" },
                            { emoji: "\uf059", label: "Q&A", slug: "q-a" },
                            { emoji: "\ue209", label: "Show and tell", slug: "show-and-tell" }
                        ]

                        Rectangle {
                            required property var modelData

                            radius: height / 2
                            color: pillMouse.containsMouse ? Qt.rgba(1, 1, 1, 0.2) : Qt.rgba(1, 1, 1, 0.1)
                            border.color: Qt.rgba(1, 1, 1, 0.25)

                            width: pill_layout.implicitWidth + 16
                            height: pill_layout.implicitHeight + 8

                            RowLayout {
                                id: pill_layout
                                anchors.fill: parent
                                Label {
                                    id: pillLabel
                                    font.family: "Font Awesome 7 Free"
                                    text: modelData.emoji
                                    color: "white"
                                }

                                Label {
                                    text: modelData.label
                                    color: "white"
                                }
                            }

                            MouseArea {
                                id: pillMouse
                                anchors.fill: parent
                                cursorShape: Qt.PointingHandCursor
                                hoverEnabled: true
                                onClicked: Qt.openUrlExternally(
                                    "https://github.com/NatLabRockies/SolTrace/discussions/categories/" + modelData.slug)
                            }
                        }


                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "GitHub"
            collapsed: false

            GridLayout {
                rows: 2
                columns: 2

                STButton {
                    text: "SolTrace"
                    text_icon: "\uf08e"

                    onClicked: {
                        Qt.openUrlExternally("https://github.com/NatLabRockies/SolTrace")
                    }
                }

                STButton {
                    text: "SolTrace Issue Tracker"
                    text_icon: "\uf08e"

                    onClicked: {
                        Qt.openUrlExternally("https://github.com/NatLabRockies/SolTrace/issues")
                    }
                }

                STButton {
                    text: "SolTrace GUI"
                    text_icon: "\uf08e"

                    onClicked: {
                        Qt.openUrlExternally("https://github.com/nicholasbl/SolTrace")
                    }
                }

                STButton {
                    text: "SolTrace GUI Issue Tracker"
                    text_icon: "\uf08e"

                    onClicked: {
                        Qt.openUrlExternally("https://github.com/nicholasbl/SolTrace/issues")
                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Build Info"
            collapsed: false

            ColumnLayout {
                spacing: 6

                Label {
                    Layout.fillWidth: true
                    text: AppData.current_version_info
                    font.bold: true
                    elide: Text.ElideRight
                }

                STClickableLabel {
                    Layout.fillWidth: true
                    text: "Copy build info to clipboard"
                    color: Material.accentColor
                    borderWidth: 0

                    onClicked: AppData.copy_build_info_to_clipboard()
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Inline Documentation"
            collapsed: false

            ColumnLayout {
                STComboBox {
                    currentIndex: App.docs.locale
                    model: ["English", "Spanish"]
                    onCurrentIndexChanged: App.docs.locale = currentIndex
                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
