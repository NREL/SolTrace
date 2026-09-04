import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Item {
    id: root
    property var blur_source

    onWidthChanged: {
        App.view.fit_panels(width, false, true)
        bottom_bar.collapsed = root.width < bottom_bar.normal_width
    }

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

    FullPanel {
        anchors.top: top_bar.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: bottom_bar.top
        anchors.margins: 10
        blur_source: root.blur_source
        enabled: App.view.full_panel.visible
        opacity: enabled
        available_width: root.width
    }

    SplitView {
        id: split
        anchors.top: top_bar.bottom
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: bottom_bar.top

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

            SplitView.minimumWidth: 400
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
            SplitView.minimumWidth: App.view.left_panel.visible && App.view.right_panel.visible ? 10 : 0
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

    BottomBar {
        id: bottom_bar
        available_width: root.width
        blur_source: root.blur_source
    }
}
