import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

import SolTrace

GridLayout {
    id: root

    property string label: ""
    property int labelTextFormat: Text.PlainText
    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property int labelWidth: 120
    property int labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    default property alias content: contentRow.data

    Layout.fillWidth: true
    columns: singleColumn ? 1 : 2
    columnSpacing: 8
    rowSpacing: 5

    STPropertyLabel {
        text: root.label
        textFormat: root.labelTextFormat
        Layout.preferredWidth: root.singleColumn ? implicitWidth : root.labelWidth
        Layout.fillWidth: root.singleColumn
        Layout.alignment: root.labelAlignment
    }

    RowLayout {
        id: contentRow

        Layout.fillWidth: true
        spacing: 8
    }
}
