import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Dialogs
import QtQuick.Layouts

import SolTrace

Flickable {
    id: root

    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property int labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    function displayPath(urlValue) {
        var text = String(urlValue)
        if (!text.length) return "Choose directory"
        if (text.startsWith("file:///")) {
            text = decodeURIComponent(text.substring(8))
            if (text.length > 1 && text[1] !== ":") {
                text = "/" + text
            }
        } else if (text.startsWith("file://")) {
            text = decodeURIComponent(text.substring(7))
        }
        return text
    }

    contentWidth: width
    contentHeight: content_column.implicitHeight
    clip: true
    boundsBehavior: Flickable.StopAtBounds

    ColumnLayout {
        id: content_column
        width: root.width

        InlineDocumentation {
            key: "analyze.export"
        }

        STPropertyPanel {
            Layout.fillWidth: true
            columns: root.singleColumn ? 1 : 2
            collapsible: true
            title: "Export Destination"

            STPropertyLabel {
                text: "Directory"
                Layout.alignment: root.labelAlignment
            }

            STButton {
                Layout.fillWidth: true
                text: root.displayPath(AppData.exporter.export_directory)
                left_text_icon: "\uf07c"
                onClicked: folderDialog.open()

                FolderDialog {
                    id: folderDialog
                    currentFolder: AppData.exporter.export_directory
                    selectedFolder: AppData.exporter.export_directory
                    onAccepted: AppData.exporter.export_directory = selectedFolder
                }
            }

            STPropertyLabel {
                text: "Result"
                Layout.alignment: root.labelAlignment
            }

            Label {
                Layout.fillWidth: true
                text: AppData.exporter.current_result_name.length
                      ? AppData.exporter.current_result_name
                      : "None"
                elide: Text.ElideRight
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true
            columns: root.singleColumn ? 1 : 2
            collapsible: true
            title: "Content"

            STPropertyLabel {
                text: "Include flux map images"
                Layout.alignment: root.labelAlignment
            }

            STSwitch {
                checked: AppData.exporter.export_flux_map_images
                text: checked ? "Enabled" : "Disabled"
                onToggled: AppData.exporter.export_flux_map_images = checked
            }

            STPropertyLabel {
                text: "Generated maps"
                Layout.alignment: root.labelAlignment
            }

            Label {
                Layout.fillWidth: true
                text: AppData.exporter.generated_flux_map_count
            }

            STPropertyLabel {
                text: "Include rays"
                Layout.alignment: root.labelAlignment
            }

            STSwitch {
                checked: AppData.exporter.export_rays
                text: checked ? "Enabled" : "Disabled"
                onToggled: AppData.exporter.export_rays = checked
            }

            STPropertyLabel {
                text: "Randomly sample rays"
                enabled: AppData.exporter.export_rays
                Layout.alignment: root.labelAlignment
            }

            STSwitch {
                enabled: AppData.exporter.export_rays
                checked: AppData.exporter.random_sample_rays
                text: checked ? "Enabled" : "Disabled"
                onToggled: AppData.exporter.random_sample_rays = checked
            }

            STPropertyLabel {
                text: "Sample size"
                visible: AppData.exporter.export_rays
                         && AppData.exporter.random_sample_rays
                Layout.alignment: root.labelAlignment
            }

            STDoubleSpinBox {
                visible: AppData.exporter.export_rays
                         && AppData.exporter.random_sample_rays
                Layout.fillWidth: true
                from: 1
                to: 1000000000
                decimals: 0
                stepSize: 1000
                suffix: "rays"
                value: AppData.exporter.random_sample_ray_count
                onValueModified: AppData.exporter.random_sample_ray_count = value
            }
        }

        STButton {
            Layout.fillWidth: true
            text: "Export"
            left_text_icon: "\uf019"
            enabled: AppData.exporter.can_export
            onClicked: AppData.exporter.export_current()
        }
    }
}
