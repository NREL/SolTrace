import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ColumnLayout {
    id: root
    spacing: 0

    Flickable {
        Layout.fillWidth: true
        Layout.fillHeight: true
        contentHeight: content.implicitHeight
        clip: true

        ColumnLayout {
            id: content
            width: parent.width
            spacing: 0

            ShadowedRectangle {
                Layout.fillWidth: true
                Layout.preferredHeight: banner.height
                Layout.topMargin: 10
                Layout.bottomMargin: 20
                radius: 8

                ColumnLayout {
                    id: banner
                    readonly property int verticalPadding: 16
                    readonly property int horizontalPadding: 20

                    width: parent.width

                    Label {
                        Layout.topMargin: banner.verticalPadding
                        Layout.leftMargin: banner.horizontalPadding
                        Layout.rightMargin: banner.horizontalPadding
                        font.pointSize: 18
                        font.family: "CMU Serif"
                        text: "Welcome to SolTrace"
                        font.bold: true
                        font.capitalization: Font.SmallCaps
                    }

                    Label {
                        Layout.fillWidth: true
                        Layout.leftMargin: banner.horizontalPadding
                        Layout.rightMargin: banner.horizontalPadding
                        Layout.bottomMargin: banner.verticalPadding
                        text: "An optical ray tracing tool for modeling concentrating solar power systems. Define your optical geometry, trace rays through the system, and analyze the resulting flux distributions."
                        wrapMode: Text.WordWrap
                    }
                }
            }

            Repeater {
                model: [
                    {
                        step: "1",
                        title: "Load Scene",
                        description: "Add a scene to the environment by loading a .stinput file. You can load multiple scenes and switch between them."
                    },
                    {
                        step: "2",
                        title: "Configure Scene",
                        description: "Set the sun shape and position, define optical properties, and arrange stages and elements that make up your concentrator geometry."
                    },
                    {
                        step: "3",
                        title: "Trace Scene",
                        description: "Configure ray tracing parameters, run the simulation, and review diagnostic output such as logs and bounding boxes."
                    },
                    {
                        step: "4",
                        title: "Analyze Results",
                        description: "Explore simulation results, inspect ray intersections, analyze flux distributions, and export data. Run multiple traces to compare results."
                    }
                ]

                delegate: ColumnLayout {
                    Layout.fillWidth: true
                    spacing: 0

                    required property var modelData

                    RowLayout {
                        Layout.fillWidth: true
                        Layout.leftMargin: 10
                        Layout.rightMargin: 10
                        spacing: 8

                        ShadowedRectangle {
                            Layout.preferredWidth: 25
                            Layout.preferredHeight: 25
                            radius: 100
                            glassColor: App.theme.glassColor

                            Label {
                                anchors.centerIn: parent
                                text: modelData.step
                            }
                        }

                        Label {
                            text: modelData.title
                            font.pointSize: 14
                            font.bold: true
                        }
                    }

                    Label {
                        Layout.fillWidth: true
                        Layout.leftMargin: 43
                        Layout.rightMargin: 10
                        Layout.topMargin: 4
                        Layout.bottomMargin: 16
                        text: modelData.description
                        wrapMode: Text.WordWrap
                    }
                }
            }

            Item { Layout.preferredHeight: 32 }
        }
    }

    Rectangle {
        Layout.fillWidth: true
        height: 1
        color: Material.dividerColor
    }

    WorkflowStepper {
        Layout.fillWidth: true
        next: "Load Scene"
        currentIndex: ViewModule.Start
    }
}
