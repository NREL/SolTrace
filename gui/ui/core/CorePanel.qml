import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Item {
    id: root
    property var blur_source

    onWidthChanged: App.view.fit_panels(width, false, true)

    // C++ -> QML (consume width from ViewModule)
    Binding {
        target: left_panel
        property: "SplitView.preferredWidth"
        value: App.view.left_panel.width
        when: left_panel.visible
    }

    Binding {
        target: right_panel
        property: "SplitView.preferredWidth"
        value: App.view.right_panel.width
        when: right_panel.visible
    }

    // QML -> C++ (write back drag-resized widths, then let fit_panels snap)
    Connections {
        target: left_panel
        function onWidthChanged() {
            if (!split.resizing) return
            if (left_panel.width === App.view.left_panel.width) return
            App.view.left_panel.width = left_panel.width
        }
    }

    Connections {
        target: right_panel
        function onWidthChanged() {
            if (!split.resizing) return
            if (right_panel.width === App.view.right_panel.width) return
            App.view.right_panel.width = right_panel.width
        }
    }

    TopBar {
        id: top_bar
        anchors.left: parent.left
        anchors.margins: 10
        anchors.right: parent.right
        anchors.top: parent.top
        blur_source: root.blur_source
        height: 42
        available_width: root.width
    }

    SettingsPanel {
        anchors.top: top_bar.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: bottom_blank.top
        anchors.margins: 10
        blur_source: root.blur_source
        enabled: App.view.settings_panel.visible
        opacity: enabled
        available_width: root.width
    }

    SplitView {
        id: split
        anchors.top: top_bar.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: bottom_blank.top
        anchors.margins: 10
        orientation: Qt.Horizontal

        handle: Rectangle {
            implicitWidth: 6
            color: SplitHandle.pressed ? "#80ffffff"
                                       : SplitHandle.hovered ? "#40ffffff"
                                                             : "transparent"
        }

        LeftPanel {
            id: left_panel

            blur_source: root.blur_source
            available_width: root.width

            SplitView.minimumWidth: 250
            SplitView.maximumWidth: split.width
            SplitView.fillHeight: true

            visible: App.view.left_panel.visible
            enabled: visible
            opacity: enabled
        }

        Item {
            id: center_space
            SplitView.fillWidth: true
            SplitView.fillHeight: true
            SplitView.minimumWidth: App.view.left_panel.visible && App.view.right_panel.visible ? 40 : 0
        }

        RightPanel {
            id: right_panel

            blur_source: root.blur_source
            available_width: root.width

            SplitView.minimumWidth: 250
            SplitView.maximumWidth: split.width
            SplitView.fillHeight: true

            visible: App.view.right_panel.visible
            enabled: visible
            opacity: enabled
        }
    }

    Item {
        id: bottom_blank

        anchors.bottom: parent.bottom
        anchors.bottomMargin: 10

        height: 42
    }

    ShadowedGlassRectangle {
        id: workflow_bar

        property bool is_open: false

        anchors.horizontalCenter: parent.horizontalCenter
        anchors.bottom: parent.bottom
        anchors.bottomMargin: 10

        blur_source: root.blur_source
        radius: 42 / 2
        glassColor: App.theme.glassColor

        width: Math.min(Math.max(parent.width
                                 - 2 * (navigation_settings_button.implicitWidth
                                        + 28), 0),
                        workflow_bar.is_open
                        ? workflow_large_pane.implicitWidth
                        : workflow_data_pane.implicitWidth)
        height: workflow_bar.is_open ? workflow_large_pane.implicitHeight : 42

        onIs_openChanged: {
            if (is_open) {
                workflow_data_pane.visible = false
                workflow_data_pane.opacity = 0
                workflow_large_pane.opacity = 0
                workflow_large_pane.visible = true
                workflow_large_fade.restart()
            } else {
                workflow_large_pane.visible = false
                workflow_large_pane.opacity = 0
                workflow_data_pane.opacity = 0
                workflow_data_pane.visible = true
                workflow_small_fade.restart()
            }
        }

        Item {
            NumberAnimation {
                id: workflow_small_fade
                target: workflow_data_pane
                property: "opacity"
                to: 1
                duration: 150
            }

            NumberAnimation {
                id: workflow_large_fade
                target: workflow_large_pane
                property: "opacity"
                to: 1
                duration: 150
            }
        }

        WorkflowPane {
            id: workflow_data_pane

            anchors.fill: parent

            opacity: 1
            visible: true

            onOpenClicked: workflow_bar.is_open = true
        }

        WorkflowPaneLarge {
            id: workflow_large_pane

            anchors.fill: parent

            opacity: 0
            visible: false

            onCloseClicked: workflow_bar.is_open = false
        }


    }

    // Rectangle {
    //     height: 32

    //     anchors.left: workflow_bar.left
    //     anchors.right: workflow_bar.right
    //     anchors.verticalCenter: workflow_bar.top

    //     anchors.verticalCenterOffset: -8

    //     border.width: 1

    //     color: open_panel_rect.containsMouse ? Material.rippleColor : "transparent"

    //     MouseArea {
    //         id: open_panel_rect

    //         anchors.fill: parent

    //         hoverEnabled: true
    //         onClicked: workflow_bar.is_open = !workflow_bar.is_open
    //     }
    // }

    STIconButton {
        id: navigation_settings_button

        anchors.right: parent.right
        anchors.rightMargin: 14
        anchors.verticalCenter: bottom_blank.verticalCenter

        text: "\uf03d"
        toolTip: "Navigation Settings"
        label.font.pointSize: 20

        onClicked: nav_settings_pop.open()

        NavigationSettings {
            id: nav_settings_pop
        }
    }
}
