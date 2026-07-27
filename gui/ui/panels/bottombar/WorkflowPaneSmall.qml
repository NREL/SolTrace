import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ShadowedGlassRectangle {
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

        spacing: 4

        Connections {
            target: AppData.file_source

            function onRowsInserted(parent, first, last) {
                root.last_db_count = AppData.file_source.rowCount()
                if (last >= first)
                    root.flash_added_data()
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

            ScriptAction { script: root.highlighted = true }
            PauseAnimation { duration: 400 }
            ScriptAction { script: root.highlighted = false }
            PauseAnimation { duration: 400 }
        }

        STIconButton {
            icon: "\uf005"
            toolTip: "Get Started"
            opacity: App.view.workflow_phase === ViewModule.Start ? 1.0 : 0.5

            onClicked: {
                App.view.simulation_content_view = false
                App.view.workflow_phase = ViewModule.Start
            }
        }

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            text: "\uf178"
            opacity: App.view.workflow_phase === ViewModule.Start
                     || App.view.workflow_phase === ViewModule.Load ? 1.0 : 0.5
        }

        STIconButton {
            icon: "\uf07c"
            toolTip: "Load Scene"
            opacity: App.view.workflow_phase === ViewModule.Load ? 1.0 : 0.5

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Load)
                    file_menu.open()
                else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Load
                }
            }

            WorkflowFileMenu {
                id: file_menu
            }
        }

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            text: "\uf178"
            opacity: App.view.workflow_phase === ViewModule.Load
                     || App.view.workflow_phase === ViewModule.Configure ? 1.0 : 0.5
        }

        STIconButton {
            icon: "\uf1de"
            toolTip: "Configure Scene"
            opacity: App.view.workflow_phase === ViewModule.Configure ? 1.0 : 0.5

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Configure)
                    data_pop.open()
                else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Configure
                }
            }
        }

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            text: "\uf178"
            opacity: App.view.workflow_phase === ViewModule.Configure
                     || App.view.workflow_phase === ViewModule.Simulate ? 1.0 : 0.5
        }

        STIconButton {
            icon: "\uf04b"
            toolTip: "Trace Scene"
            opacity: App.view.workflow_phase === ViewModule.Simulate ? 1.0 : 0.5

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Simulate)
                    data_pop.open()
                else {
                    App.view.simulation_content_view = false
                    App.view.workflow_phase = ViewModule.Simulate
                }
            }
        }

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            text: "\uf178"
            opacity: App.view.workflow_phase === ViewModule.Simulate
                     || App.view.workflow_phase === ViewModule.Analyze ? 1.0 : 0.5
        }

        STIconButton {
            icon: "\uf080"
            toolTip: "Analyze Results"
            opacity: App.view.workflow_phase === ViewModule.Analyze ? 1.0 : 0.5

            onClicked: {
                if (App.view.workflow_phase === ViewModule.Analyze)
                    results_pop.open()
                else {
                    App.view.simulation_content_view = true
                    App.view.workflow_phase = ViewModule.Analyze
                }
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
