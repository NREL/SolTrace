import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root

    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth
    clip: true

    readonly property var thirdPartyLicenses: loadThirdPartyLicenses()

    function loadThirdPartyLicenses() {
        const request = new XMLHttpRequest()
        request.open("GET", "qrc:/licenses/third_party_licenses.json", false)
        request.send()
        if (request.status !== 200 && request.status !== 0) {
            return {
                "licenses": [],
                "warnings": [qsTr("Embedded third-party license manifest could not be loaded.")]
            }
        }

        try {
            return JSON.parse(request.responseText)
        } catch (error) {
            return {
                "licenses": [],
                "warnings": [qsTr("Embedded third-party license manifest could not be parsed.")]
            }
        }
    }

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: qsTr("Licenses")
        }

        Label {
            Layout.fillWidth: true
            wrapMode: Text.Wrap
            font.pixelSize: App.theme._normalSize
            color: App.theme.fontColor
            text: qsTr("SolTrace uses the following third-party assets. All skybox environment maps are licensed under CC0 1.0 Universal (Public Domain) and require no attribution.")
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: qsTr("Skybox: Clear Sky")
            collapsed: true

            ColumnLayout {
                spacing: 6

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("Autumn Field Pure Sky")
                    font.bold: true
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    textFormat: Text.MarkdownText
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("**Source:** [Poly Haven](https://polyhaven.com/a/autumn_field_puresky)  \n**Authors:** Jarod Guest (Sky Edits), Sergej Majboroda (Original)  \n**License:** [CC0 1.0 Universal (Public Domain)](https://polyhaven.com/license)")
                    onLinkActivated: (link) => Qt.openUrlExternally(link)

                    MouseArea {
                        anchors.fill: parent
                        acceptedButtons: Qt.NoButton
                        cursorShape: parent.hoveredLink ? Qt.PointingHandCursor : Qt.ArrowCursor
                    }
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    opacity: 0.7
                    text: qsTr("Ground recolored to desert sand. Decorative only — does not affect ray-trace results.")
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: qsTr("Skybox: Partly Cloudy")
            collapsed: true

            ColumnLayout {
                spacing: 6

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("Kloofendal 48d Partly Cloudy Pure Sky")
                    font.bold: true
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    textFormat: Text.MarkdownText
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("**Source:** [Poly Haven](https://polyhaven.com/a/kloofendal_48d_partly_cloudy_puresky)  \n**Authors:** Greg Zaal (Original), Jarod Guest (Sky Edits)  \n**License:** [CC0 1.0 Universal (Public Domain)](https://polyhaven.com/license)")
                    onLinkActivated: (link) => Qt.openUrlExternally(link)

                    MouseArea {
                        anchors.fill: parent
                        acceptedButtons: Qt.NoButton
                        cursorShape: parent.hoveredLink ? Qt.PointingHandCursor : Qt.ArrowCursor
                    }
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    opacity: 0.7
                    text: qsTr("Synthetic desert sand ground added to lower hemisphere. Decorative only — does not affect ray-trace results.")
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            title: qsTr("Skybox: Sunset")
            collapsed: true

            ColumnLayout {
                spacing: 6

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("Industrial Sunset 02 Pure Sky")
                    font.bold: true
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    textFormat: Text.MarkdownText
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    text: qsTr("**Source:** [Poly Haven](https://polyhaven.com/a/industrial_sunset_02_puresky)  \n**Authors:** Jarod Guest (Sky Edits), Sergej Majboroda (Original)  \n**License:** [CC0 1.0 Universal (Public Domain)](https://polyhaven.com/license)")
                    onLinkActivated: (link) => Qt.openUrlExternally(link)

                    MouseArea {
                        anchors.fill: parent
                        acceptedButtons: Qt.NoButton
                        cursorShape: parent.hoveredLink ? Qt.PointingHandCursor : Qt.ArrowCursor
                    }
                }

                Label {
                    Layout.fillWidth: true
                    wrapMode: Text.Wrap
                    font.pixelSize: App.theme._normalSize
                    color: App.theme.fontColor
                    opacity: 0.7
                    text: qsTr("Ground recolored to desert sand. Decorative only — does not affect ray-trace results.")
                }
            }
        }

        Repeater {
            model: root.thirdPartyLicenses.warnings || []

            Label {
                required property string modelData

                Layout.fillWidth: true
                wrapMode: Text.Wrap
                font.pixelSize: App.theme._normalSize
                color: Material.color(Material.Yellow)
                text: modelData
            }
        }

        Repeater {
            model: root.thirdPartyLicenses.licenses || []

            STPropertyPanel {
                id: licensePanel

                required property var modelData

                Layout.fillWidth: true
                title: licensePanel.modelData.component + ": " + licensePanel.modelData.name
                collapsed: true
                columns: 1

                ColumnLayout {
                    Layout.fillWidth: true
                    spacing: 6

                    Label {
                        Layout.fillWidth: true
                        wrapMode: Text.WrapAnywhere
                        font.pixelSize: App.theme._normalSize
                        color: App.theme.fontColor
                        opacity: 0.7
                        text: licensePanel.modelData.path
                    }

                    ScrollView {
                        id: lic_scroll
                        Layout.fillWidth: true
                        Layout.preferredHeight: Math.min(360, Math.max(120, implicitHeight))

                        contentWidth: availableWidth

                        TextArea {
                            width: lic_scroll.availableWidth
                            readOnly: true
                            selectByMouse: false
                            wrapMode: TextEdit.Wrap
                            textFormat: TextEdit.PlainText
                            text: licensePanel.modelData.text
                            color: App.theme.fontColor
                            font.family: "monospace"
                            font.pixelSize: App.theme._normalSize
                        }
                    }


                }
            }
        }

        Item { Layout.fillHeight: true }
    }
}
