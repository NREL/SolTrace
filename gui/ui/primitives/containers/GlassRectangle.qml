import QtQuick
import QtQuick.Controls
import QtQuick.Effects
import SolTrace

Rectangle {
    id: root

    property var blur_source
    property int blurMax: 32
    property color glassColor: App.theme.glassColor

    clip: true
    color: "transparent"

    function getPositionRelativeToWindow() {
        try {
            var pos = root.mapToItem(blur_source, 0, 0)
            return Qt.rect(pos.x, pos.y, root.width, root.height)
        } catch (e) {
            console.error("Error getting position:", e)
        }
        return Qt.rect(0, 0, 0, 0)
    }


    Rectangle {
        id: mask
        anchors.fill: parent
        layer.enabled: true
        radius: root.radius
    }


    ShaderEffectSource {
        id: effect_source
        anchors.fill: root
        sourceItem: root.blur_source
        visible: false
        hideSource: false
        live: !root.isDestroying

        sourceRect: {
            // We track changes to the size of the root window
            root.width; root.height; root.x; root.y

            if (root.blur_source) {
                root.blur_source.width; root.blur_source.height
            }

            let p = root.parent
            while (p) { p.x; p.y; p.width; p.height; p = p.parent }

            if (!root.blur_source) return Qt.rect(0, 0, 0, 0)
            const pos = root.mapToItem(root.blur_source, 0, 0)
            return Qt.rect(pos.x, pos.y, root.width, root.height)
        }
    }


    MultiEffect {
        id: multiEffect
        anchors.fill: effect_source
        source: effect_source
        maskEnabled: true
        maskSource: mask
        autoPaddingEnabled: false
        blurEnabled: true
        blurMax: root.blurMax
        blurMultiplier: 1
        blur: 1
    }


    Rectangle {
        id: background
        color: root.glassColor
        border.width: 1
        border.color: Qt.rgba(0.9, 0.9, 0.9, 0.15)
        anchors.fill: parent
        radius: root.radius
    }
}
