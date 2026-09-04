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

        InlineDocumentation {
            key: "settings.overview"
            alwaysVisible: true
            showTitle: true
            Layout.fillWidth: true
        }
        InlineDocumentation {
            key: "settings.background"
            alwaysVisible: true
            showTitle: true
            Layout.fillWidth: true
        }
        InlineDocumentation {
            key: "settings.methodology"
            alwaysVisible: true
            showTitle: true
            Layout.fillWidth: true
        }
    }
}
