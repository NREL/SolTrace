import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    Label {
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8

        text: "Explore simulation results, inspect ray intersections, analyze flux distributions, and export data."
        wrapMode: Text.WordWrap
    }


    STComboBar {
        id: bar

        currentIndex: Math.min(App.view.analyze_section, 3)
        onCurrentIndexChanged: App.view.analyze_section = currentIndex

        Layout.fillWidth: true

        collapseLabels: App.view.left_panel.size === SplitPanelData.Small

        iconModel: ["\uf03a", "\ue4bc", "\uf06d", "\uf019"]
        model: ["Results", "Intersections", "Flux", "Export"]
    }

    Rectangle {
        Layout.fillWidth: true
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    StackLayout {
        currentIndex: bar.currentIndex
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.margins: 8

        ResultListPane {}

        AnalyzeIntersections {

        }

        AnalyzeFlux {

        }

        AnalyzeExport {

        }
    }

    Rectangle {
        Layout.fillWidth: true
        height: 1
        color: Material.dividerColor
    }

    WorkflowStepper {
        previous: "Run Tracer"
        currentIndex: ViewModule.Analyze
    }
}
