import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Controls.Material

import SolTrace

ShadowedGlassRectangle {
    id: root

    required property int available_width

    signal show_script_area()

    readonly property bool show_logos: available_width > 900
    radius: height / 2
    glassColor: App.theme.glassColor

    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 10
        anchors.rightMargin: 10
        spacing: 0

        TopBarLeftPane {
            id: top_bar_left_pane

            Layout.fillHeight: true
            Layout.fillWidth: true

            available_width: root.available_width
            show_logos: root.show_logos
        }

        Item {
            Layout.fillWidth: true
        }

        TopBarRightPane {
            Layout.fillHeight: true

            blur_source: root.blur_source
            available_width: root.available_width
        }
    }
}
