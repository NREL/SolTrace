import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root
    Layout.fillWidth: true

    InlineDocumentation {
        key: "workflow.trace"
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8
        showTitle: false
    }

    Rectangle {
        Layout.fillWidth: true
        Layout.preferredHeight: 1
        color: Material.dividerColor
    }

    TracePanel {
        Layout.fillWidth: true
        Layout.fillHeight: true
    }

    Rectangle {
        Layout.fillWidth: true
        height: 1
        color: Material.dividerColor
    }

    WorkflowStepper {
        previous: "Configure Scene"
        next: "Analyze Results"
        currentIndex: ViewModule.Simulate
    }
}
