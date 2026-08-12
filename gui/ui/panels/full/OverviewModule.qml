import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: "Overview"
        }

        InlineDocumentation {
            key: "settings.overview"
            alwaysVisible: true
            showTitle: false
            Layout.fillWidth: true
        }
    }
}
