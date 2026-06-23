import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

RowLayout {
    id: root

    required property int available_width
    required property bool show_logos

    readonly property bool analyzing: App.view.workflow_phase === 3
    readonly property string current_context_title: analyzing ? "Analyze" : "Configure"
    readonly property string current_context_name: analyzing
                                                   ? (App.simulation.current_simulation_result_name.length ?
                                                          App.simulation.current_simulation_result_name
                                                        : "None")
                                                   : (App.file_source.current_database ?
                                                          App.file_source.current_database.name
                                                        : "None")

    spacing: 0

    STIconButton {
        id: leftpanel_open
        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitHeight
        //Layout.leftMargin: 20
        //Layout.rightMargin: 10

        text: "\uf0c9"
        toolTip: (App.view.left_panel.visible ? "Close": "Open") + " Left Panel"

        label.font.pointSize: 20
        onClicked: {
            if (App.view.settings_panel.visible) {
                App.view.settings_panel.visible = false
                return
            }

            App.view.left_panel.visible = !App.view.left_panel.visible && !App.view.settings_panel.visible
            App.view.fit_panels(root.available_width)
        }
    }

    Item {
        id: soltrace_logo
        visible: root.show_logos

        Layout.alignment: Qt.AlignVCenter

        Layout.rightMargin: 10

        implicitWidth: logo_content.implicitWidth
        implicitHeight: logo_content.implicitHeight

        Rectangle {
            anchors.fill: logo_content
            anchors.margins: -5
            radius: 12

            color: logo_mouse_area.containsMouse ? Material.rippleColor :
                                                   "transparent"

            Behavior on color {
                ColorAnimation {
                    duration: 100
                }
            }
        }

        Column {
            id: logo_content
            spacing: -2

            RowLayout {
                id: logo_row
                spacing: 4

                Item { Layout.preferredWidth: 2 }

                Image {
                    id: logo_icon
                    source: "qrc:/assets/images/logo.svg"
                    Layout.preferredHeight: 27
                    Layout.preferredWidth: Layout.preferredHeight
                    Layout.alignment: Qt.AlignBottom
                    mipmap: true
                    fillMode: Image.PreserveAspectFit

                    layer.enabled: true
                    layer.effect: MultiEffect {
                        colorization: 1.0
                        colorizationColor: App.theme.fontColor
                    }
                }

                Label {
                    id: logo_text
                    Layout.fillWidth: true
                    Layout.alignment: Qt.AlignBottom
                    font.pointSize: 18
                    font.family: "CMU Serif"
                    text: "SolTrace"
                    font.bold: true
                    font.capitalization: Font.SmallCaps
                }

                STClickableLabel {
                    id: prerel_logo

                    visible: AppData.is_prerelease

                    //anchors.bottom: parent.bottom
                    //anchors.left: parent.left

                    text: "\uf071"
                    font.pointSize: 24
                    font.family: "Font Awesome 7 Free"

                    style: Label.Outline
                    styleColor: "black"

                    color: Material.color(Material.Yellow)
                }

                Item { Layout.preferredWidth: 2 }
            }

            Rectangle {
                color: App.theme.fontColor
                width: logo_row.width
                height: 1
            }
        }

        MouseArea {
            id: logo_mouse_area
            anchors.fill: logo_content
            hoverEnabled: true
            onClicked: version_pop.open()
        }

        STPopup {
            id: version_pop

            contentWidth: 280

            ColumnLayout {
                anchors.fill: parent
                spacing: 6

                Label {
                    Layout.fillWidth: true

                    text: "SolTrace"
                    font.bold: true
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                }

                Label {
                    Layout.fillWidth: true

                    color: Material.color(Material.Yellow)

                    text: "This version of SolTrace is intended for testing and evaluation. Features, file formats, and simulation behavior may change before the final release. Verify important results with a stable release before using them for production work."
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                }

                Label {
                    Layout.fillWidth: true

                    text: AppData.current_build_info
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                }
            }
        }
    }

    RowLayout {
        id: context_row

        Layout.fillWidth: true
        Layout.leftMargin: 8
        Layout.rightMargin: 12
        spacing: 8

        Label {
            text: root.current_context_title
            font.bold: true
            opacity: 0.65
            elide: Text.ElideRight
        }

        Rectangle {
            Layout.fillHeight: true
            Layout.topMargin: 10
            Layout.bottomMargin: 10
            Layout.preferredWidth: 1
            color: Material.dividerColor
        }

        Label {
            Layout.fillWidth: true
            Layout.maximumWidth: 320
            text: root.current_context_name
            elide: Text.ElideMiddle
        }
    }

}
