import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

Item {
    id: root

    implicitWidth: mode_row.implicitWidth + mode_row.anchors.leftMargin
                   + mode_row.anchors.rightMargin
    implicitHeight: mode_row.implicitHeight

    property int last_db_count: AppData.file_source.rowCount()
    property bool highlighted: false
    function flash_added_data() {
        flash_highlight_animation.restart()
    }

    signal openClicked()

    RowLayout {
        id: mode_row

        anchors.fill: parent
        anchors.leftMargin: 20
        anchors.rightMargin: 20

        spacing: 8

        Connections {
            target: AppData.file_source

            function onRowsInserted(parent, first, last) {
                root.last_db_count = AppData.file_source.rowCount()
                if (last >= first) {
                    root.flash_added_data()
                }
            }

            function onRowsRemoved(parent, first, last) {
                root.last_db_count = AppData.file_source.rowCount()
                flash_highlight_animation.stop()
                root.highlighted = false
            }
        }

        SequentialAnimation {
            id: flash_highlight_animation

            loops: 3

            ScriptAction {
                script: root.highlighted = true
            }

            PauseAnimation {
                duration: 400
            }

            ScriptAction {
                script: root.highlighted = false
            }

            PauseAnimation {
                duration: 400
            }
        }

        Label {
            text: "Workflow"

            //font.family: "CMU Serif"

            font.bold: true

            opacity: 0.65

            Rectangle {
                id: wf_highlight

                bottomLeftRadius: 42 / 2
                topLeftRadius: 42 / 2
                anchors.fill: parent
                anchors.leftMargin: -20
                anchors.rightMargin: -8
                anchors.topMargin: -12
                anchors.bottomMargin: -12

                color: Material.rippleColor

                opacity: wf_mouse.containsMouse

                Behavior on opacity {
                    NumberAnimation {
                        duration: 150
                    }
                }
            }

            MouseArea {
                id: wf_mouse
                anchors.fill: parent
                anchors.margins: -5

                hoverEnabled: true

                onClicked: root.openClicked()
            }
        }

        Rectangle {
            color: Material.dividerColor

            width: 1
            Layout.fillHeight: true
        }

        STClickableLabel {
            id: data_label

            Layout.fillHeight: true
            Layout.leftMargin: 4
            Layout.rightMargin: 4

            borderWidth: 0

            property bool is_active: App.view.workflow_phase === 0

            text: "1. Data"
            font.pointSize: 16
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter
            opacity: is_active ? 1.0 : .50

            onClicked: {
                if (App.view.workflow_phase === 0) {
                    file_menu.open()
                } else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = 0
                }
            }

            WorkflowFileMenu {
                id: file_menu
            }
        }

        Label {
            id: file_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 4
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf101"

            property bool highlight: data_label.is_active
                                     || current_scene_label.is_active

            opacity: highlight ? 1.0 : .50
        }

        STClickableLabel {
            id: current_scene_label

            property bool is_active: App.view.workflow_phase === 1

            Layout.fillHeight: true

            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            text: "2. Configure"
            elide: Label.ElideMiddle

            //font.bold: is_active
            font.pointSize: 16
            opacity: is_active ? 1.0 : .50

            onClicked: {
                if (App.view.workflow_phase === 1) {
                    data_pop.open()
                } else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = 1
                }
            }

            Behavior on opacity {
                NumberAnimation {
                    duration: 150
                }
            }
        }

        STClickableLabel {
            id: swap_data_indicator
            text: "\uf0d7"

            font.family: "Font Awesome 7 Free"

            visible: App.view.workflow_phase === 1

            onClicked: data_pop.open()
        }

        Label {
            id: configure_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 4
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf101"

            property bool highlight: current_scene_label.is_active || simulate_label.is_active

            opacity: highlight ? 1.0 : .50

        }

        STClickableLabel {
            id: simulate_label

            property bool is_active: App.view.workflow_phase === 2

            Layout.fillHeight: true

            text: "3. Trace"
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            //font.bold: is_active
            font.pointSize: 16
            opacity: is_active ? 1.0 : .50

            onClicked: {
                App.view.simulation_content_view = false
                App.view.workflow_phase = 2
            }

            Behavior on opacity {
                NumberAnimation {
                    duration: 150
                }
            }
        }

        STClickableLabel {
            id: swap_simulation_data_indicator
            text: "\uf0d7"

            font.family: "Font Awesome 7 Free"

            visible: App.view.workflow_phase === 2

            onClicked: data_pop.open()
        }

        Label {
            id: analyze_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 4
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf101"

            property bool highlight: analyze_label.is_active || simulate_label.is_active

            opacity: highlight ? 1.0 : .50
        }

        STClickableLabel {
            id: analyze_label

            property bool is_active: App.view.workflow_phase === 3

            Layout.fillHeight: true

            font.pointSize: 16

            text: "4. Analyze"
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            //font.bold: is_active
            opacity: is_active ? 1.0 : .50

            onClicked: {
                if (App.view.workflow_phase === 3) {
                    results_pop.open()
                } else {
                    App.view.simulation_content_view = true
                    App.view.workflow_phase = 3
                }
            }

            Behavior on opacity {
                NumberAnimation {
                    duration: 150
                }
            }
        }

        STClickableLabel {
            id: swap_results_indicator
            text: "\uf0d7"

            font.family: "Font Awesome 7 Free"

            visible: App.view.workflow_phase === 3

            onClicked: results_pop.open()
        }
    }

    Item {
        anchors.fill: parent

        ResultsPopup {
            id: results_pop

            //y: root.popup_above ? -height - 10 : root.height + 10
            width: root.width
            height: Overlay.overlay.height * 0.66
        }


        DataPopup {
            id: data_pop

            //y: root.popup_above ? -height - 10 : root.height + 10
            width: root.width
            height: Overlay.overlay.height * 0.66
        }
    }
}
