import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Dialogs
import QtQuick.Layouts

import SolTrace

Flickable {
    id: root

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

        STFormPanel {
            collapsible: true
            title: "Export Destination"

            STFormRow {
                label: "Directory"

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
            }

            STFormRow {
                label: "Result"

                Label {
                    Layout.fillWidth: true
                    text: AppData.exporter.current_result_name.length
                          ? AppData.exporter.current_result_name
                          : "None"
                    elide: Text.ElideRight
                }
            }
        }

        STFormPanel {
            collapsible: true
            title: "Content"

            STFormRow {
                label: "Generated maps"

                Label {
                    Layout.fillWidth: true
                    text: AppData.exporter.generated_flux_map_count
                }
            }

            STSwitch {
                checked: AppData.exporter.export_flux_map_images
                text: "Include flux map images"
                onToggled: AppData.exporter.export_flux_map_images = checked
            }

            STSwitch {
                checked: AppData.exporter.export_flux_map_meshes
                text: "Include flux map meshes"
                onToggled: AppData.exporter.export_flux_map_meshes = checked
            }

            STSwitch {
                checked: AppData.exporter.export_rays
                text: "Include rays"
                onToggled: AppData.exporter.export_rays = checked
            }

            STSwitch {
                enabled: AppData.exporter.export_rays
                checked: AppData.exporter.random_sample_rays
                text: "Randomly sample rays"
                onToggled: AppData.exporter.random_sample_rays = checked
            }

            STFormRow {
                label: "Sample size"
                visible: AppData.exporter.export_rays
                         && AppData.exporter.random_sample_rays

                STDoubleSpinBox {
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

            STSwitch {
                checked: AppData.exporter.export_scene_copy
                text: "Include copy of scene"
                onToggled: AppData.exporter.export_scene_copy = checked
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
