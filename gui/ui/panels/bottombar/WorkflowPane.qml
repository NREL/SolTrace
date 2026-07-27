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

    required property var blur_source
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

        spacing: 4

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

            Layout.rightMargin: 12

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

            Layout.preferredWidth: 1
            Layout.fillHeight: true
        }

        ShadowedRectangle {
            Layout.preferredWidth: 25
            Layout.preferredHeight: 25
            Layout.leftMargin: 8

            radius: 100
            blur_source: root.blur_source
            glassColor: App.theme.glassColor

            Label {
                text: "1"
                anchors.centerIn: parent
            }
        }

        STClickableLabel {
            id: start_label

            Layout.fillHeight: true
            Layout.rightMargin: 4

            borderWidth: 0

            property bool is_active: App.view.workflow_phase === ViewModule.Start

            text: App.view.workflow_phase === ViewModule.Start ? "Get Started" : "Start"
            font.pointSize: 16
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter
            opacity: is_active ? 1.0 : .50

            onClicked: {
                App.view.simulation_content_view = false
                App.view.workflow_phase = ViewModule.Start
            }
        }

        Label {
            id: start_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 8
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf178"

            property bool highlight: start_label.is_active || data_label.is_active

            opacity: highlight ? 1.0 : .50
        }

        ShadowedRectangle {
            Layout.preferredWidth: 25
            Layout.preferredHeight: 25
            Layout.leftMargin: 4

            radius: 100
            blur_source: root.blur_source
            glassColor: App.theme.glassColor

            Label {
                text: "2"
                anchors.centerIn: parent
            }
        }

        STClickableLabel {
            id: data_label

            Layout.fillHeight: true
            Layout.rightMargin: 4

            borderWidth: 0

            property bool is_active: App.view.workflow_phase === ViewModule.Load

            text: App.view.workflow_phase === ViewModule.Load ? "Load Scene" : "Load"
            font.pointSize: 16
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter
            opacity: is_active ? 1.0 : .50

            onClicked: {
                App.view.simulation_content_view = false
                App.view.workflow_phase = ViewModule.Load
            }
        }

        STClickableLabel {
            id: file_menu_indicator
            text: "\uf0d7"

            font.family: "Font Awesome 7 Free"

            visible: App.view.workflow_phase === ViewModule.Load

            onClicked: file_menu.open()

            WorkflowFileMenu {
                id: file_menu
            }
        }

        Label {
            id: file_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 8
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf178"

            property bool highlight: data_label.is_active
                                     || current_scene_label.is_active

            opacity: highlight ? 1.0 : .50
        }

        ShadowedRectangle {
            Layout.preferredWidth: 25
            Layout.preferredHeight: 25
            Layout.leftMargin: 4

            radius: 100
            blur_source: root.blur_source
            glassColor: App.theme.glassColor

            Label {
                text: "3"
                anchors.centerIn: parent
            }
        }

        STClickableLabel {
            id: current_scene_label

            property bool is_active: App.view.workflow_phase === ViewModule.Configure

            Layout.fillHeight: true

            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            text: App.view.workflow_phase === ViewModule.Configure ? "Configure Scene" : "Configure"
            elide: Label.ElideMiddle

            //font.bold: is_active
            font.pointSize: 16
            opacity: is_active ? 1.0 : .50

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Configure) {
                    data_pop.open()
                } else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Configure
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

            visible: App.view.workflow_phase === ViewModule.Configure

            onClicked: data_pop.open()
        }

        Label {
            id: configure_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 8
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf178"

            property bool highlight: current_scene_label.is_active || simulate_label.is_active

            opacity: highlight ? 1.0 : .50

        }

        ShadowedRectangle {
            Layout.preferredWidth: 25
            Layout.preferredHeight: 25
            Layout.leftMargin: 4

            radius: 100
            blur_source: root.blur_source
            glassColor: App.theme.glassColor

            Label {
                text: "4"
                anchors.centerIn: parent
            }
        }

        STClickableLabel {
            id: simulate_label

            property bool is_active: App.view.workflow_phase === ViewModule.Simulate

            Layout.fillHeight: true

            text: App.view.workflow_phase === ViewModule.Simulate ? "Trace Scene" : "Trace"
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            //font.bold: is_active
            font.pointSize: 16
            opacity: is_active ? 1.0 : .50

            onClicked: {
                App.view.simulation_content_view = false
                App.view.workflow_phase = ViewModule.Simulate
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

            visible: App.view.workflow_phase === ViewModule.Simulate

            onClicked: data_pop.open()
        }

        Label {
            id: analyze_separator

            Layout.leftMargin: 4
            Layout.rightMargin: 8
            Layout.alignment: Qt.AlignVCenter

            font.family: "Font Awesome 7 Free"
            text: "\uf178"

            property bool highlight: analyze_label.is_active || simulate_label.is_active

            opacity: highlight ? 1.0 : .50
        }

        ShadowedRectangle {
            Layout.preferredWidth: 25
            Layout.preferredHeight: 25
            Layout.leftMargin: 4

            radius: 100
            blur_source: root.blur_source
            glassColor: App.theme.glassColor

            Label {
                text: "5"
                anchors.centerIn: parent
            }
        }

        STClickableLabel {
            id: analyze_label

            property bool is_active: App.view.workflow_phase === ViewModule.Analyze

            Layout.fillHeight: true

            font.pointSize: 16

            text: App.view.workflow_phase === ViewModule.Analyze ? "Analyze Results" : "Analyze"
            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            //font.bold: is_active
            opacity: is_active ? 1.0 : .50

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Analyze) {
                    results_pop.open()
                } else {
                    App.view.simulation_content_view = true
                    App.view.workflow_phase = ViewModule.Analyze
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

            visible: App.view.workflow_phase === ViewModule.Analyze

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
