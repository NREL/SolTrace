import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material
import SolTrace

ColumnLayout {
    spacing: 8

    InlineDocumentation {
        key: "workflow.load"
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8
        showTitle: false
    }

    SceneListPane {
        Layout.fillWidth: true
        Layout.fillHeight: true
    }

    WorkflowStepper {
        previous: "Get Started"
        next: "Configure Scene"
        currentIndex: ViewModule.Load
    }
}
