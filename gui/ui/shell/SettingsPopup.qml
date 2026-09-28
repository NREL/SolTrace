import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STPopup {
    id: root

    contentWidth: Math.min(760, Math.max(360, availableWidth - 40))
    contentHeight: 560
    modal: false
    closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

    ColumnLayout {
        anchors.fill: parent
        spacing: 10

        RowLayout {
            Layout.fillWidth: true

            STIconButton {
                visible: !settingsView.wideMode && settingsView.editing
                icon: "\uf060"
                toolTip: "Back to Settings"

                onClicked: settingsView.goBack()
            }

            Label {
                Layout.fillWidth: true

                text: "Settings"
                font.bold: true
                font.pointSize: App.theme.labelSize + 2
                wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            }

            STIconButton {
                icon: "\uf00d"

                onClicked: root.close()
            }
        }

        AdaptiveEditor {
            id: settingsView

            Layout.fillWidth: true
            Layout.fillHeight: true

            wideThreshold: 500
            listWidth: 200
            currentIndex: App.view.settings_section
            onCurrentIndexChanged: {
                App.view.settings_section = currentIndex
            }

            model: ListModel {
                ListElement { name: "General"; icon: "\uf013" }
                ListElement { name: "Theme"; icon: "\uf53f" }
                ListElement { name: "Scene"; icon: "\uf1b2" }
                ListElement { name: "Language"; icon: "\uf0ac" }
            }

            listFooter: ColumnLayout {
                spacing: 8

                STDangerousButton {
                    Layout.fillWidth: true
                    text: "Reset All Settings"
                    left_text_icon: "\uf2ed"
                    onClicked: resetSettingsDialog.open()
                }

                STDialog {
                    id: resetSettingsDialog
                    title: "Reset All Settings"
                    modal: true

                    ColumnLayout {
                        anchors.fill: parent

                        Label {
                            Layout.fillWidth: true
                            text: "This will reset all saved settings to their defaults."
                            wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                        }
                    }

                    footer: STDialogButtonBox {
                        STButton {
                            text: "Cancel"
                            DialogButtonBox.buttonRole: DialogButtonBox.RejectRole
                        }

                        STDangerousButton {
                            text: "Reset Settings"
                            DialogButtonBox.buttonRole: DialogButtonBox.DestructiveRole
                            onClicked: {
                                AppData.reset_all_settings()
                                resetSettingsDialog.close()
                            }
                        }
                    }
                }
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
                currentIndex: settingsView.currentIndex

                GeneralSettings {}

                ThemeSettings {}

                SceneSettings {}

                LanguageSettings {}
            }
        }
    }
}
