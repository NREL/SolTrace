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
    //color: "transparent"
    //border.width: 1
    //border.color: Theme.lineColor
    //border.color: Material.dividerColor
    opacity: enabled ? 1.0 : 0.55
    implicitHeight: form_core.implicitHeight
    implicitWidth: form_core.implicitWidth

    clip: true

    default property alias contentChildren: layout.children

    ColumnLayout {
        id: form_core

        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.bottom: parent.bottom

        spacing: 0

        Item {
            implicitHeight: header_row.implicitHeight + 12
            Layout.preferredHeight: implicitHeight
            //implicitWidth: header_row.implicitWidth + 12

            //topLeftRadius: root.radius
            //topRightRadius: root.radius

            Layout.fillWidth: true
            //Layout.bottomMargin: root.collapsed ? 8 : 0

            visible: root.title.length > 0

            //color: Material.dividerColor

            RowLayout {
                id: header_row

                anchors.fill: parent
                anchors.margins: 6

                Label {
                    font.family: "Font Awesome 7 Free"
                    visible: root.collapsible
                    text: root.collapsed ? "\uf105" : "\uf107"


                }

                CheckBox {
                    visible: root.checkable
                    id: check_box
                }

                Label {
                    text: root.title
                    Layout.fillWidth: true

                    font.pointSize: App.theme.propertyPanelHeaderSize

                }

                // Rectangle {
                //     Layout.fillWidth: true
                //     Layout.preferredHeight: 1
                //     color: Material.dividerColor
                // }

            }

            MouseArea {
                anchors.fill: parent
                onClicked: root.collapsed = !root.collapsed
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
            Layout.fillHeight: true

            Layout.margins: 6

            columns: root.columns

            visible: !root.collapsed
        }
    }
}
