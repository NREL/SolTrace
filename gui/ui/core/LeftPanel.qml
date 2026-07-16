import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ShadowedGlassRectangle {
    id: root

    property int available_width: 0

    radius: 10

    MouseArea {
        anchors.fill: parent
        acceptedButtons: Qt.AllButtons
        hoverEnabled: true
        preventStealing: true

        onPressed: (mouse) => mouse.accepted = true
        onReleased: (mouse) => mouse.accepted = true
        onClicked: (mouse) => mouse.accepted = true
        onDoubleClicked: (mouse) => mouse.accepted = true
        onPositionChanged: (mouse) => mouse.accepted = true
        onWheel: (wheel) => wheel.accepted = true
    }

    RowLayout {
        id: module_info_row
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 10
        height: 34

        Item {
            Layout.fillWidth: true
            visible: !App.view.left_panel.is_small()
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.minimumWidth: 0
            spacing: App.view.left_panel.is_small() ? 6 : 10

            Item {
                Layout.preferredWidth: 5
                visible: App.view.left_panel.is_small()
            }

            ShadowedRectangle {
                Layout.preferredWidth: 25
                Layout.preferredHeight: 25
                Layout.leftMargin: 4

                radius: 100
                blur_source: root.blur_source
                glassColor: App.theme.glassColor

                visible: App.view.workflow_phase > 0

                Label {
                    text: App.view.workflow_phase
                    anchors.centerIn: parent
                }
            }

            Label {
                text: ["Get Started", "Load Scene", "Configure Scene", "Trace Scene", "Analyze Results"][
                          Math.min(App.view.workflow_phase, 4)]
                font.pointSize: 16
                font.bold: true
                font.family: "CMU Serif"
            }

            /*
            Repeater {
                visible: false
                model: ["Data", "Configure", "Simulate", "Analyze"]

                RowLayout {
                    id: labelRow
                    required property int index
                    required property string modelData
                    property var icons: ["\uf0ad", "\uf085", "\uf201"]

                    property bool is_active : App.view.workflow_phase === index

                    spacing: 4

                    STClickableLabel {
                        text: parent.icons[parent.index]
                        font.family: "Font Awesome 7 Free"
                        opacity: is_active ? 1 : 0.5
                        borderWidth: 0
                        font.pointSize: App.view.left_panel.is_small() ? 16 : 12

                        onClicked: {
                            App.view.workflow_phase = parent.index
                            App.view.simulation_content_view = parent.index === 3
                        }

                        STToolTip {
                            visible: parent.containsMouse && App.view.left_panel.size == SplitPanelData.Small
                            text: labelRow.modelData
                        }
                    }

                    STClickableLabel {
                        text: parent.modelData
                        borderWidth: 0
                        font.pointSize: 16
                        font.bold: true
                        font.family: "CMU Serif"
                        font.underline: is_active && !App.view.left_panel.is_small()
                        opacity: is_active ? 1 : 0.5
                        visible: {
                            if (is_active && App.view.left_panel.width >= 220) {
                                return true
                            }

                            if (App.view.left_panel.width >= 450) {
                                return true
                            }

                            return false

                        }

                        onClicked: {
                            App.view.workflow_phase = parent.index
                            App.view.simulation_content_view = parent.index === 3
                        }
                    }

                    Label {
                        Layout.leftMargin: 4
                        Layout.rightMargin: 4
                        font.family: "Font Awesome 7 Free"
                        text: "\uf101"
                        visible: parent.index < 3
                    }
                }
            }*/

            Item {
                Layout.fillWidth: true
                visible: !App.view.left_panel.is_small()
            }
        }

        Item {
            Layout.fillWidth: true
            visible: !App.view.left_panel.is_small()
        }

        PanelButtons {
            id: panel_buttons
            Layout.alignment: Qt.AlignRight | Qt.AlignVCenter
            Layout.preferredWidth: implicitWidth
            target: App.view.left_panel
            otherTarget: App.view.right_panel
            available_width: root.available_width
            is_right_panel: false
        }
    }

    StackLayout {
        id: module_stack
        currentIndex: App.view.workflow_phase
        onCurrentIndexChanged: App.view.workflow_phase = currentIndex

        anchors.top: module_info_row.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 10
        anchors.leftMargin: 16
        anchors.rightMargin: 16

        StartModule {}
        LoadModule {}
        ConfigureModule {}
        SimulateModule {}
        AnalyzeModule {}
    }
}
