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

                Label {
                    text: Math.max(1, App.view.workflow_phase)
                    anchors.centerIn: parent
                }
            }

            Label {
                text: ["Load Scene", "Configure Scene", "Trace Scene", "Analyze Results"][
                          Math.max(0, Math.min(App.view.workflow_phase - 1, 3))]
                font.pointSize: 16
                font.bold: true
                font.family: "CMU Serif"
            }

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
        currentIndex: Math.max(0, App.view.workflow_phase - 1)
        onCurrentIndexChanged: App.view.workflow_phase = currentIndex + 1

        anchors.top: module_info_row.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 10
        anchors.leftMargin: 16
        anchors.rightMargin: 16

        LoadModule {}
        ConfigureModule {}
        SimulateModule {}
        AnalyzeModule {}
    }
}
