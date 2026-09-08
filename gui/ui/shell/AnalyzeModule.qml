import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    InlineDocumentation {
        key: "workflow.analyze"
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8
        showTitle: false
    }


    STComboBar {
        id: bar

        currentIndex: Math.min(App.view.analyze_section, 3)
        onCurrentIndexChanged: App.view.analyze_section = currentIndex

        Layout.fillWidth: true

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
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    WorkflowStepper {
        previous: "Trace Scene"
        currentIndex: ViewModule.Analyze
    }
}
