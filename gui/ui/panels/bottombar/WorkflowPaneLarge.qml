import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

Item {
    id: root

    implicitWidth: workflow_layout.implicitWidth + 32
    implicitHeight: workflow_layout.implicitHeight + 18

    readonly property string active_name: App.file_source.current_database ?
                                             App.file_source.current_database.name : "None"
    readonly property string active_result_name: App.simulation.current_simulation_result_name

    signal closeClicked()

    ColumnLayout {
        id: workflow_layout

        anchors.fill: parent
        anchors.leftMargin: 16
        anchors.rightMargin: 16
        anchors.topMargin: 8
        anchors.bottomMargin: 10
        spacing: 8

        Rectangle {
            Layout.fillWidth: true

            implicitHeight: wf_label.implicitHeight + 12

            radius: 6

            opacity: 0.65

            color: wf_mouse.containsMouse ? Material.rippleColor : "transparent"

            Label {
                id: wf_label
                text: "Workflow"
                font.bold: true

                anchors.fill: parent

                verticalAlignment: Qt.AlignVCenter

                Label {
                    anchors.verticalCenter: parent.verticalCenter
                    anchors.right: parent.right
                    font.family: "Font Awesome 7 Free"
                    text: "\uf0d7"
                }
            }

            MouseArea {
                id: wf_mouse
                anchors.fill: parent

                hoverEnabled: true

                onClicked: root.closeClicked()
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 10

            WorkflowLargeItem {
                icon: "\uf56f"
                title: "1. Load"
                value: "Load"
                active: App.view.workflow_phase === ViewModule.Load

                description: "Load scene"

                onActivated: {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Load
                }

                onValueClicked: file_menu.open()

                WorkflowFileMenu {
                    id: file_menu
                }
            }

            Label {
                Layout.alignment: Qt.AlignVCenter
                font.family: "Font Awesome 7 Free"
                text: "\uf101"
                opacity: App.view.workflow_phase === ViewModule.Load
                         || App.view.workflow_phase === ViewModule.Configure ? 1.0 : 0.5
            }

            WorkflowLargeItem {
                icon: "\uf7d9"
                title: "2. Configure"
                value: root.active_name
                active: App.view.workflow_phase === ViewModule.Configure

                description: "Arrange and define the trace scenario"

                onActivated: {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Configure
                }

                onValueClicked: data_pop.open()
            }

            Label {
                Layout.alignment: Qt.AlignVCenter
                font.family: "Font Awesome 7 Free"
                text: "\uf101"
                opacity: App.view.workflow_phase === ViewModule.Configure
                         || App.view.workflow_phase === ViewModule.Simulate ? 1.0 : 0.5
            }

            WorkflowLargeItem {
                icon: "\uf04b"
                title: "3. Trace"
                value: root.active_name
                active: App.view.workflow_phase === ViewModule.Simulate

                description: "Execute a ray trace of a scene"

                onActivated: {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Simulate
                }

                onValueClicked: data_pop.open()
            }

            Label {
                Layout.alignment: Qt.AlignVCenter
                font.family: "Font Awesome 7 Free"
                text: "\uf101"
                opacity: App.view.workflow_phase === ViewModule.Simulate
                         || App.view.workflow_phase === ViewModule.Analyze ? 1.0 : 0.5
            }

            WorkflowLargeItem {
                icon: "\uf1fe"
                title: "4. Analyze"
                value: root.active_result_name
                active: App.view.workflow_phase === ViewModule.Analyze

                description: "Analyze traced rays"

                onActivated: {
                    App.view.simulation_content_view = true
                    App.view.workflow_phase = ViewModule.Analyze
                }

                onValueClicked: results_pop.open()
            }
        }
    }

    Item {
        anchors.fill: parent

        ResultsPopup {
            id: results_pop

            width: root.width
            height: Overlay.overlay.height * 0.66
        }

        DataPopup {
            id: data_pop

            width: root.width
            height: Overlay.overlay.height * 0.66
        }
    }
}
