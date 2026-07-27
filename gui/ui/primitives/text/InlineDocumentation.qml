import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root
    required property string key

    property string titlePrefix: ""
    property string title: ""
    property string body: ""
    property bool showTitle: true
    property bool alwaysVisible: false
    property var blocks: []
    readonly property int docsVersion: App.docs.version
    readonly property int docsLocale: App.docs.locale

    visible: root.alwaysVisible || App.view.inline_docs

    function svgColor(color) {
        var value = String(color)
        return value.length === 9 && value[0] === "#" ? "#" + value.slice(3) : value
    }

    function svgImageSource(svg) {
        var color = svgColor(App.theme.fontColor)
        var patchedSvg = svg
            .replace(/fill="#000000"/gi, "fill=\"" + color + "\"")
            .replace(/fill="#000"/gi, "fill=\"" + color + "\"")
            .replace(/stroke="#000000"/gi, "stroke=\"" + color + "\"")
            .replace(/stroke="#000"/gi, "stroke=\"" + color + "\"")
        return "data:image/svg+xml;utf8," + encodeURIComponent(patchedSvg)
    }

    function refreshBlocks() {
        // make sure this refreshes when locale or version changes
        App.docs.locale
        App.docs.version

        if (root.body) {
            root.blocks = [{ "type": "text", "content": root.body }]
            return
        }

        try {
            root.blocks = JSON.parse(App.docs.blocks(root.key))
        } catch (error) {
            root.blocks = [{
                "type": "text",
                "content": App.docs.get(root.key)
            }]
        }
    }

    Component.onCompleted: refreshBlocks()
    onKeyChanged: refreshBlocks()
    onBodyChanged: refreshBlocks()
    onDocsVersionChanged: refreshBlocks()
    onDocsLocaleChanged: refreshBlocks()

    SubHeader {
        id: title
        visible: root.showTitle
        text: {
            App.docs.locale
            App.docs.version
            return root.title ? root.title : App.docs.get(root.key, "title")
        }

        Layout.fillWidth: true
    }

    Repeater {
        model: root.blocks

        delegate: Loader {
            id: blockLoader

            property var block: modelData

            Layout.fillWidth: true
            sourceComponent: block.type === "math" ? mathBlock : textBlock

            Component {
                id: textBlock

                Label {
                    width: blockLoader.width
                    text: blockLoader.block.content
                    textFormat: Label.MarkdownText
                    wrapMode: Text.WordWrap
                    font.pointSize: App.theme.normalSize
                    color: App.theme.fontColor
                }
            }

            Component {
                id: mathBlock

                Item {
                    width: blockLoader.width
                    implicitHeight: mathImage.implicitHeight

                    Image {
                        id: mathImage
                        anchors.horizontalCenter: parent.horizontalCenter
                        source: root.svgImageSource(blockLoader.block.svg)
                    }
                }
            }
        }
    }
}
