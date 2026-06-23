import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts
import SolTrace

Flow {
    id: root
    property var model: null
    property Component delegate: null
    property Component preview: null
    spacing: 8

    // Propagate model data
    function injectProperties(item, index) {
        const data = root.model.get(index)
        for (const key of Object.keys(data)) {
            if (item.hasOwnProperty(key))
                item[key] = data[key]
        }
    }

    Repeater {
        model: root.model
        delegate: Rectangle {
            id: wrapper
            required property int index

            width: 300
            height: cardContent.implicitHeight + 24
            color: App.theme.glassColorA(0.1)
            border.color: Qt.alpha("black", .25)
            radius: 10

            ColumnLayout {
                id: cardContent
                width: parent.width - 32
                anchors.centerIn: parent
                spacing: 5

                Loader {
                    Layout.fillWidth: true
                    sourceComponent: root.delegate
                    onLoaded: root.injectProperties(item, wrapper.index)
                }
            }

            STIconButton {
                anchors.top: parent.top
                anchors.right: parent.right
                anchors.margins: 8
                text: "\uf08e"
                onClicked: popup.open()
                visible: root.preview !== null
            }

            Popup {
                id: popup
                width: 400
                modal: true
                anchors.centerIn: Overlay.overlay
                Overlay.modal: Rectangle {
                    color: Qt.rgba(0, 0, 0, 0.25)
                }
                padding: 12
                closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

                Loader {
                    anchors.left: parent.left
                    anchors.right: parent.right
                    sourceComponent: root.preview
                    onLoaded: root.injectProperties(item, wrapper.index)
                }
            }
        }
    }
}
