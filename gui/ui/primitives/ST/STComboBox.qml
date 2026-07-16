pragma ComponentBehavior: Bound

import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

ComboBox {
    id: root
    property int lastIndex: currentIndex < 0 ? 0 : currentIndex

    Layout.preferredWidth: Math.max(implicitWidth + 50, 200)
    Material.foreground: App.theme.fontColor
    font.pointSize: App.theme.labelSize

    onCurrentIndexChanged: {
        if (currentIndex < 0) currentIndex = lastIndex
        else lastIndex = currentIndex
    }

    background: WellRectangle {
        implicitWidth: 80
        implicitHeight: 32
        radius: height / 2
    }

    indicator: Label {
        font.family: "Font Awesome 7 Free"
        text: "\uf078"

        x: root.width - width - root.rightPadding
        y: root.topPadding + (root.availableHeight - height) / 2
        width: 12
        height:12
    }

    popup: Popup {
        //y: root.height - 1
        width: root.width
        height: Math.min(contentItem.implicitHeight, root.Window.height - topMargin - bottomMargin)
        padding: 2

        //onOpened: console.log("popup opened")
        //onClosed: console.log("popup closed")

        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

        contentItem: ListView {
            clip: true
            implicitHeight: contentHeight
            model: root.popup.visible ? root.delegateModel : null
            currentIndex: root.highlightedIndex

            ScrollBar.vertical: STScrollBar { }

            delegate: STItemDelegate {
                id: delegate
                width: ListView.view.width
                text: root.textAt(index)

                required property var model
                required property int index


                background: Rectangle {
                    implicitHeight: 24
                    implicitWidth: 100
                    opacity: enabled ? 1 : 0.3
                    color: parent.down
                           ? Material.rippleColor : "transparent"
                }

                contentItem: Label {
                    text: delegate.text
                    font: delegate.font
                }
            }
        }

        background: Rectangle {
            color: Qt.alpha(Material.backgroundColor, .90)
            border.color: Material.frameColor
            radius: 10
        }
    }
}
