import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Material

import SolTrace

RowLayout {
    id: root
    property string title

    Layout.columnSpan: parent.columns === 2 ? 2 : 1

    Label {
        text: root.title

        visible: text.length > 0
    }

    Rectangle {
        Layout.fillWidth: true
        color: Material.dividerColor
        height: 1
    }

}
