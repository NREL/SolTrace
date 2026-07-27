import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts
import QtQuick.Dialogs

import SolTrace

ScrollView {
    Layout.fillHeight: true
    Layout.fillWidth: true

    id: root

    property var module: AppData.script

    contentWidth: availableWidth

    function load_script(url) {
        const requestUrl = String(url).startsWith(":/")
                         ? "qrc" + String(url)
                         : url

        const request = new XMLHttpRequest()
        request.open("GET", requestUrl)
        request.onreadystatechange = function() {
            if (request.readyState === XMLHttpRequest.DONE) {
                if (request.status === 0 || request.status === 200) {
                    root.module.code = request.responseText
                } else {
                    console.warn("Unable to load script", requestUrl, request.status)
                    root.module.notify_error("Could not load the script.")
                }
            }
        }
        request.send()
    }

    function pathToFileUrl(path) {
        return "file://" + path
    }

    ListModel {
        id: output_model
    }

    Connections {
        target: AppData.script

        function onLogged(type, msg) {
            output_model.append({
                type: type,
                message: msg,
            })
        }
    }

    STPopup {
        id: api_popup
        parent: Overlay.overlay
        modal: true
        focus: true
        enable_shadow: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

        property string apiText: ""
        readonly property real popupMargin: 24

        onOpened: apiText = root.module.api_markdown()

        width: Math.max(320, Math.min(parent.width - popupMargin * 2, 920))
        height: Math.max(320, Math.min(parent.height - popupMargin * 2, 720))
        x: Math.round((parent.width - width) / 2)
        y: Math.round((parent.height - height) / 2)

        contentItem: ColumnLayout {
            spacing: 8

            RowLayout {
                Layout.fillWidth: true

                Label {
                    Layout.fillWidth: true

                    text: "Script API"
                    font.bold: true
                    elide: Label.ElideRight
                }

                STIconButton {
                    icon: "\uf00d"
                    toolTip: "Close"
                    onClicked: api_popup.close()
                }
            }

            ScrollView {
                Layout.fillWidth: true
                Layout.fillHeight: true

                clip: true

                TextArea {
                    text: api_popup.apiText
                    readOnly: true
                    selectByMouse: true
                    wrapMode: TextEdit.Wrap
                    textFormat: TextEdit.MarkdownText

                    background: Item {}
                }
            }
        }
    }

    ColumnLayout {
        width: root.availableWidth

        InlineDocumentation {
            key: "script.interface"
            Layout.fillWidth: true
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Source"
            collapsible: true

            RowLayout {
                uniformCellSizes: true

                Layout.fillWidth: true
                Layout.columnSpan: 2

                STButton {
                    Layout.fillWidth: true

                    text: "New"
                }

                STButton {
                    Layout.fillWidth: true

                    text: "Load"

                    onClicked: file_menu.open()

                    STMenu {
                        id: file_menu
                        MenuItem {
                            text: "Open"
                            onClicked: openFileDialog.open()
                        }

                        STMenu {
                            id: recents_menu
                            title: qsTr("Builtins")
                            enabled: recent_instantiator.count > 0

                            Instantiator {
                                id: recent_instantiator
                                model: root.module.builtin_scripts

                                delegate: MenuItem {
                                    implicitWidth: 140
                                    text: modelData.replace(":/assets/scripts/", "")
                                    onTriggered: root.load_script(modelData)
                                }
                                onObjectAdded: (index, object) => recents_menu.insertItem(index, object)
                                onObjectRemoved: (index, object) => recents_menu.removeItem(object)
                            }
                        }
                    }

                    FileDialog {
                        id: openFileDialog
                        onAccepted: {
                            var str_file = String(selectedFile)

                            currentFolder = str_file.substring(0, str_file.lastIndexOf("/"))
                            root.load_script(str_file)
                        }
                    }
                }

                STButton {
                    Layout.fillWidth: true

                    text: "Edit"

                    onClicked: script_editor.open()

                    ScriptEditor {
                        id: script_editor
                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: "Parameters"
            collapsible: true

            Label {
                Layout.fillWidth: true
                Layout.columnSpan: 2

                Layout.row: 0

                elide: Label.ElideRight

                property bool has_title: root.module.title.length

                text: has_title ? root.module.title : "Available script parameters will appear here."
                font.bold: has_title
            }

            Label {
                visible: root.module.description.length
                Layout.fillWidth: true
                Layout.columnSpan: 2

                Layout.row: 1

                wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                text: root.module.description
            }

            Repeater {
                id: label_repeater
                model: root.module.properties

                delegate: STPropertyLabel {
                    required property string name
                    required property int index

                    Layout.column: 0
                    Layout.row: index + 2
                    Layout.fillWidth: true
                    text: name
                }
            }

            Repeater {
                model: root.module.properties

                delegate: STTextField {
                    required property string value
                    required property int index
                    required property var model

                    Layout.column: 1
                    Layout.row: index + 2
                    Layout.fillWidth: true

                    text: value

                    onAccepted: model.value = text
                    onTextEdited: model.value = text
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true
            title: "Execute"

            RowLayout {
                Layout.fillWidth: true
                Layout.columnSpan: 2


                STTextField {
                    Layout.fillWidth: true
                    text: root.module.working_directory

                    onAccepted: root.module.working_directory = text
                    onTextEdited: root.module.working_directory = text
                }

                // something weird with this icon button. doesnt seem to have the right bounds
                STIconButton {
                    icon: "\uf07c"
                    toolTip: "Select Working Directory"
                    onClicked: workingDirectoryDialog.open()


                    FolderDialog {
                        id: workingDirectoryDialog
                        currentFolder: root.pathToFileUrl(root.module.working_directory)

                        onAccepted: {
                            const value = String(selectedFolder)
                            root.module.working_directory =
                                    decodeURIComponent(value.replace(/^file:\/\//, ""))
                        }
                    }
                }
            }

            STButton {
                text: "Run"
                Layout.columnSpan: 2
                Layout.fillWidth: true

                enabled: root.module.valid

                onClicked: {
                    output_model.clear()
                    root.module.run()
                }
            }

        }


        STPropertyPanel {
            Layout.fillWidth: true

            title: "Output"
            collapsible: true

            ListView {
                Layout.fillWidth: true
                Layout.columnSpan: 2

                Layout.preferredHeight: 150

                model: output_model
                clip: true

                ScrollBar.vertical: STScrollBar { }

                delegate: Label {
                    required property string message
                    required property int type

                    width: ListView.view.width
                    height: implicitHeight

                    text: message

                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere

                    color: {
                        switch (type) {
                        case 0: Material.foreground; break;
                        case 1: Material.color(Material.Yellow); break;
                        case 2: Material.color(Material.Red); break;
                        }
                    }
                }

                onCountChanged: positionViewAtEnd()
            }
        }

        STIconButton {
            Layout.alignment: Qt.AlignRight
            Layout.rightMargin: 6

            icon: "\uf02d"
            label: "Script API"
            toolTip: "Show Script API"
            onClicked: api_popup.open()
        }
    }
}
