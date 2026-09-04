// TODO make this work.

import QtQuick
import QtQuick.Layouts

Item {
    id: root

    default property alias content: container.data

    property int animationDuration: 300

    readonly property Item contentItem:
        container.children.length > 0 ? container.children[0] : null

    implicitWidth: contentItem ? contentItem.implicitWidth : 0
    implicitHeight: contentItem ? contentItem.implicitHeight : 0

    // Give the Layout sensible preferred dimensions.
    Layout.preferredWidth: implicitWidth
    Layout.preferredHeight: implicitHeight

    clip: true

    Behavior on x {
        NumberAnimation {
            duration: root.animationDuration
            easing.type: Easing.InOutCubic
        }
    }

    Behavior on y {
        NumberAnimation {
            duration: root.animationDuration
            easing.type: Easing.InOutCubic
        }
    }

    Item {
        id: container
        anchors.fill: parent
    }
}