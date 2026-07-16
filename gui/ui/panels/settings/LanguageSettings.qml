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
            text: "Language"
        }

        ColumnLayout {
            STComboBox {
                currentIndex: App.docs.locale
                model: ["English", "Spanish"]
                onCurrentIndexChanged: App.docs.locale = currentIndex
            }
        }
    }
}

