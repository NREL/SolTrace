import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Material
import SolTrace

Rectangle {
    id: root
    property string title

    property bool collapsible: true
    property bool checkable: false

    property bool collapsed: false
    property alias checked: check_box.checked

    property int columns: 2

    Behavior on implicitHeight {
        NumberAnimation {
            duration: 50
        }
    }

    radius: 14
    color: Qt.rgba(1, 1, 1, 0.03)
    border.width: 1
    //border.color: Theme.lineColor
    border.color: Material.dividerColor
    opacity: enabled ? 1.0 : 0.55
    implicitHeight: form_core.implicitHeight + 20
    implicitWidth: form_core.implicitWidth + 20

    clip: true

    // onHeightChanged: console.log(height)
    // onImplicitHeightChanged: console.log(implicitHeight)

    default property alias contentChildren: layout.children

    ColumnLayout {
        id: form_core

        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 10

        spacing: 10

        RowLayout {
            visible: root.title.length > 0
            Layout.bottomMargin: root.collapsed ? 8 : 0

            Label {
                visible: root.collapsible
                text: root.collapsed ? "▼" : "▶"

                MouseArea {
                    anchors.fill: parent
                    onClicked: root.collapsed = !root.collapsed
                }
            }

            CheckBox {
                visible: root.checkable
                id: check_box
            }

            Label {
                text: root.title
                font.bold: true

                font.pointSize: App.theme.propertyPanelHeaderSize
                
            }
        }

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 1
            color: Material.dividerColor

            visible: !root.collapsed && root.title.length > 0
        }

        GridLayout {
            id: layout
            Layout.fillWidth: true

            columns: root.columns

            visible: !root.collapsed
        }
    }
}
