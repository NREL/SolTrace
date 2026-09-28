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

    readonly property var tintOptions: [
        { label: "Default", value: App.theme.glassTintDefault },
        { label: "Red", value: App.theme.glassTintRed },
        { label: "Pink", value: App.theme.glassTintPink },
        { label: "Purple", value: App.theme.glassTintPurple },
        { label: "Indigo", value: App.theme.glassTintIndigo },
        { label: "Blue", value: App.theme.glassTintBlue },
        { label: "Cyan", value: App.theme.glassTintCyan },
        { label: "Teal", value: App.theme.glassTintTeal },
        { label: "Green", value: App.theme.glassTintGreen },
        { label: "Lime", value: App.theme.glassTintLime },
        { label: "Yellow", value: App.theme.glassTintYellow },
        { label: "Amber", value: App.theme.glassTintAmber },
        { label: "Orange", value: App.theme.glassTintOrange },
        { label: "Brown", value: App.theme.glassTintBrown },
        { label: "Grey", value: App.theme.glassTintGrey },
    ]

    function reset() {
        App.theme.glassTint = App.theme.glassTintDefault
        App.theme.glassClarity = App.theme.glassClarityBalanced
        App.theme.fontSize = App.theme.fontSizeDefault
    }

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Theme"
        }

        STPropertyPanel {
            Layout.fillWidth: true
            title: "Glass Tint"
            collapsed: false

            Flow {
                Layout.fillWidth: true
                Layout.columnSpan: 2
                spacing: 8

                Repeater {
                    model: root.tintOptions

                    Rectangle {
                        id: tintPill
                        required property var modelData

                        readonly property bool selected: App.theme.glassTint === modelData.value

                        width: tintRow.implicitWidth + 18
                        height: Math.max(30, tintRow.implicitHeight + 8)
                        radius: height / 2
                        color: selected ? Qt.alpha(Material.accentColor, 0.18)
                                        : App.theme.glassColorA(0.08)
                        border.width: selected ? 2 : 1
                        border.color: selected ? Material.accentColor
                                               : Material.dividerColor

                        RowLayout {
                            id: tintRow
                            anchors.centerIn: parent
                            spacing: 7

                            Rectangle {
                                Layout.preferredWidth: 14
                                Layout.preferredHeight: 14
                                radius: 7
                                color: App.theme.tintColor(tintPill.modelData.value)
                                border.color: Qt.alpha(App.theme.fontColor, 0.3)
                            }

                            Label {
                                text: tintPill.modelData.label
                            }
                        }

                        MouseArea {
                            anchors.fill: parent
                            cursorShape: Qt.PointingHandCursor
                            hoverEnabled: true
                            onClicked: App.theme.glassTint = tintPill.modelData.value
                        }
                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true
            title: "Glass Translucency"
            collapsed: false

            STComboBar {
                Layout.fillWidth: true
                Layout.columnSpan: 2
                model: ["Clear", "Balanced", "Solid"]
                currentIndex: App.theme.glassClarity
                onCurrentIndexChanged: App.theme.glassClarity = currentIndex
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true
            title: "Font Size"
            collapsed: false

            STComboBar {
                Layout.fillWidth: true
                Layout.columnSpan: 2
                model: ["Small", "Default", "Large"]
                currentIndex: App.theme.fontSize
                onCurrentIndexChanged: App.theme.fontSize = currentIndex
            }
        }

        Item { Layout.fillHeight: true }
    }
}
