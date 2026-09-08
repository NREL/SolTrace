import QtQuick
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var module: App.layout.instance_edit
    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property int labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    ElementAssignmentEditor {
        Layout.fillWidth: true
        module: root.module
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment
    }

    ElementPlacementEditor {
        Layout.fillWidth: true
        module: root.module
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment
    }

    ElementVisualizationEditor {
        Layout.fillWidth: true
        module: root.module
        singleColumn: root.singleColumn
    }
}
