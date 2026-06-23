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

        console.log("Loading script from:", requestUrl)
        const request = new XMLHttpRequest()
        request.open("GET", requestUrl)
        request.onreadystatechange = function() {
            if (request.readyState === XMLHttpRequest.DONE) {
                if (request.status === 0 || request.status === 200) {
                    console.log("Script ready, installing.", request.responseText)
                    root.module.code = request.responseText
                } else {
                    console.warn("Unable to load script", requestUrl, request.status)
                    root.module.notify_error("Could not load the script.")
                }
            }
        }
        request.send()
    }

    ColumnLayout {
        width: root.availableWidth

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

            title: "Script Parameters"
            collapsible: true

            Label {
                visible: root.module.title.length
                Layout.fillWidth: true
                Layout.columnSpan: 2

                Layout.row: 0

                elide: Label.ElideRight

                text: root.module.title
                font.bold: true
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

            STButton {
                text: "Run"
                Layout.columnSpan: 2
                Layout.fillWidth: true

                enabled: root.module.valid

                onClicked: root.module.run()
            }


        }

    }
}
