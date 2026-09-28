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

    function getTitle() {
        switch (App.view.workflow_phase) {
            case ViewModule.Load: return "Scene"
            case ViewModule.Configure: return "Scene"
            case ViewModule.Simulate: return "Scene"
            case ViewModule.Analyze: return "Results"
            default: return "Other"
        }
    }

    readonly property bool analyzing: App.view.workflow_phase === ViewModule.Analyze
    readonly property bool configuring: App.view.workflow_phase === ViewModule.Configure
    readonly property string current_context_title: getTitle()
    readonly property string current_context_name: analyzing
                                                   ? (App.simulation.current_simulation_result_name.length ?
                                                          App.simulation.current_simulation_result_name
                                                        : "None")
                                                   : (App.file_source.current_database ?
                                                          App.file_source.current_database.name
                                                        : "None")

    spacing: 0

    FileController {
        id: file_controller
    }

    STIconButton {
        id: leftpanel_open
        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitHeight
        //Layout.leftMargin: 20
        //Layout.rightMargin: 10

        icon: "\uf0c9"
        toolTip: (App.view.left_panel.visible ? "Close": "Open") + " Left Panel"

        iconSize: 20
        onClicked: App.view.toggle_left_panel(root.available_width);
    }

    Item {
        id: soltrace_logo
        visible: root.show_logos

        Layout.alignment: Qt.AlignVCenter

        Layout.leftMargin: 8
        Layout.rightMargin: 15

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

        RowLayout {
            id: logo_content
            spacing: 4

            Image {
                id: logo_icon
                source: "qrc:/assets/images/logo.png"
                Layout.preferredHeight: 22
                Layout.preferredWidth: Layout.preferredHeight
                Layout.alignment: Qt.AlignCenter
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
                Layout.alignment: Qt.AlignCenter
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

                text: " \uf071"
                font.family: "Font Awesome 7 Free"

                style: Label.Outline
                styleColor: "black"

                font.pointSize: App.theme.labelSize

                color: Material.color(Material.Yellow)
            }
        }

        MouseArea {
            id: logo_mouse_area
            anchors.fill: logo_content
            hoverEnabled: true
            onClicked: info_pop.open()
        }

        STPopup {
            id: info_pop

            contentWidth: Math.min(760, Math.max(360, root.available_width - 40))
            contentHeight: 560
            modal: false
            closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

            ColumnLayout {
                anchors.fill: parent
                spacing: 10

                RowLayout {
                    Layout.fillWidth: true

                    Label {
                        Layout.fillWidth: true

                        text: "SolTrace"
                        font.bold: true
                        font.pointSize: App.theme.labelSize + 2
                        wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                    }

                    STIconButton {
                        icon: "\uf00d"

                        onClicked: info_pop.close()
                    }
                }

                Label {
                    Layout.fillWidth: true

                    color: Material.color(Material.Yellow)
                    visible: AppData.is_prerelease

                    text: "This version of SolTrace is intended for testing and evaluation. Features, file formats, and simulation behavior may change before the final release. Verify important results with a stable release before using them for production work."
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                }

                AdaptiveEditor {
                    id: infoView

                    Layout.fillWidth: true
                    Layout.fillHeight: true

                    wideThreshold: 520
                    listWidth: 180
                    currentIndex: App.view.info_section
                    onCurrentIndexChanged: {
                        App.view.info_section = currentIndex
                    }

                    model: ListModel {
                        ListElement { name: "Get Started"; icon: "\uf005" }
                        ListElement { name: "What's New"; icon: "\uf587" }
                        ListElement { name: "Overview"; icon: "\ue0bb" }
                        ListElement { name: "Community"; icon: "\uf500" }
                        ListElement { name: "Licenses"; icon: "\ue447" }
                        ListElement { name: "Build & Logs"; icon: "\uf7d9" }
                    }

                    listDelegate: ItemDelegate {
                        text: itemModel ? itemModel.name : ""
                        highlighted: isCurrent
                        width: parent ? parent.width : implicitWidth

                        contentItem: RowLayout {
                            spacing: 8
                            Label {
                                text: itemModel ? itemModel.icon : ""
                                font.family: "Font Awesome 7 Free"
                                font.pointSize: 14

                                Layout.preferredWidth: 32
                            }
                            Label {
                                text: itemModel ? itemModel.name : ""
                                Layout.fillWidth: true
                            }
                        }

                        background: Rectangle {
                            implicitHeight: 36
                            implicitWidth: 100
                            opacity: enabled ? 1 : 0.3
                            color: parent.down ? Material.rippleColor
                                 : parent.highlighted ? Qt.rgba(Material.accentColor.r,
                                                                 Material.accentColor.g,
                                                                 Material.accentColor.b, 0.12)
                                 : "transparent"
                            radius: 14
                        }
                    }

                    detailView: StackLayout {
                        currentIndex: infoView.currentIndex

                        GetStartedModule {}

                        FeatureModule {}

                        OverviewModule {}

                        CommunityModule {}

                        LicenseModule {}

                        DiagnosticsModule {}
                    }
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
            Layout.maximumWidth: 320
            text: root.current_context_name
            elide: Text.ElideMiddle
        }

        Rectangle {
            Layout.fillHeight: true
            Layout.topMargin: 10
            Layout.bottomMargin: 10
            Layout.preferredWidth: 1
            color: Material.dividerColor
        }

        STIconButton {
            visible: !root.analyzing
            enabled: !App.file_source.is_loading
            icon: "\uf56f"
            toolTip: "Load Scene"
            onClicked: file_menu.open()

            WorkflowFileMenu {
                id: file_menu
                showNewScene: false
            }
        }

        STIconButton {
            visible: !root.analyzing && !!App.file_source.current_database
            enabled: !App.file_source.is_loading
            icon: "\uf0c7"
            toolTip: "Save Scene"
            onClicked: file_controller.save_current()
        }

        Item {
            Layout.fillWidth: true
        }
    }

}
