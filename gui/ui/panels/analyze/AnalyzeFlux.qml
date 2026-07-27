import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Flickable {
    id: root
    property var left_panel_size: App.view.left_panel.size
    property var flux_module : AppData.flux

    function formatCoordinate(value) {
        if (!isFinite(value)) {
            return "--"
        }

        return Number(value).toFixed(3)
    }

    function formatCentroid(centroid) {
        if (!centroid) {
            return "--"
        }

        return "(" + formatCoordinate(centroid.x) + ", "
                   + formatCoordinate(centroid.y) + ", "
                   + formatCoordinate(centroid.z) + ")"
    }

    contentWidth: width
    contentHeight: content_column.implicitHeight
    clip: true
    boundsBehavior: Flickable.StopAtBounds

    ColumnLayout {
        id: content_column
        width: root.width

        InlineDocumentation {
            key: "analyze.flux"
        }

        STPropertyPanel {
            Layout.fillWidth: true

            collapsible: true
            title: "Current Flux Map"

            visible: AppData.flux.current_image.length

            Image {
                Layout.columnSpan: 2
                Layout.fillWidth: true
                Layout.preferredHeight: width

                source: map_selector.currentIndex == 0 ?
                            AppData.flux.current_image
                          :
                            AppData.flux.current_image + "_point_map"
                fillMode: Image.PreserveAspectFit

                mipmap: true

                STIconButton {
                    icon: "\uf019"

                    anchors.bottom: parent.bottom
                    anchors.right: parent.right
                    anchors.margins: 10
                }
            }

            STComboBar {
                id: map_selector
                Layout.columnSpan: 2
                Layout.fillWidth: true
                collapseLabels: AppData.view.left_panel.size === SplitPanelData.Small
                iconModel: ["\uf00a", "\uf141"]
                model: ["Bins", "Points"]
            }

            STSwitch {
                Layout.columnSpan: 2
                text: "Show Whole Scene"
                checked: AppData.flux.show_other_geometry
                onToggled: AppData.flux.show_other_geometry = checked
            }

            STPropertyLabel {
                text: "Plotted Power"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.plotted_power
                font.bold: true
            }

            STPropertyLabel {
                text: "Peak Flux"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.peak_flux
                font.bold: true
            }

            STPropertyLabel {
                text: "Min Flux"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.min_flux
                font.bold: true
            }

            STPropertyLabel {
                text: "Average Flux"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.average_flux
                font.bold: true
            }

            STPropertyLabel {
                text: "Sigma Flux"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.sigma_flux
                font.bold: true
            }

            STPropertyLabel {
                text: "Uniformity"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.uniformity
                font.bold: true
            }

            STPropertyLabel {
                text: "Peak Flux Uncert"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.peak_flux_uncertainty
                font.bold: true
            }

            STPropertyLabel {
                text: "Average Flux Uncert"
            }

            Label {
                Layout.fillWidth: true
                text: root.flux_module.current_flux_stats.average_flux_uncertainty
                font.bold: true
            }

            STPropertyLabel {
                text: "Centroid"
            }

            Label {
                property vector3d cent: root.flux_module.current_flux_stats.centroid
                Layout.fillWidth: true
                text: root.formatCentroid(cent)
                font.bold: true
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            collapsible: true
            title: "Compute Flux Map"

            STPropertyLabel {
                text: "Element"
            }

            STButton {
                Layout.fillWidth: true

                text: AppData.flux.current_entity_name.length ?
                          AppData.flux.current_entity_name : "Select Element"
                left_text_icon: "\uf05b"

                onClicked: entity_pop.open()

                SelectItemPopup {
                    id: entity_pop
                    source_model: AppData.flux.entity_model

                    onSelectedEntity: (entity) => {
                        AppData.flux.select_entity(entity)
                    }
                }
            }

            STButton {
                Layout.fillWidth: true
                Layout.columnSpan: 2

                text: "Compute Map"
                left_text_icon: "\uf0da"

                onClicked: {
                    AppData.flux.start_generate()
                }
            }

            STPropertySeparator {
                title: "Computed Maps"
                visible: existing_maps.count > 0
            }

            ListView {
                id: existing_maps
                visible: count > 0
                Layout.columnSpan: 2
                Layout.fillHeight: true
                Layout.fillWidth: true
                Layout.preferredHeight: 200
                clip: true

                model: AppData.flux.computed_maps_model

                ScrollBar.vertical: STScrollBar { }

                delegate: STItemDelegate {
                    id: delegate
                    required property string name
                    required property var entity
                    text: delegate.name
                    //highlighted: isCurrent

                    onClicked: {
                        AppData.flux.select_entity(delegate.entity)
                    }
                }
            }
        }

        STPropertyPanel {
            Layout.fillWidth: true

            // Disable for now until this is a working feature
            visible: false

            collapsible: true
            title: "Flux Volume"

            STPropertyLabel {
                text: "Grid Resolution"
            }

            STSpinBox {
                id: resolution_spin
                Layout.fillWidth: true

                value: 512
                from: 64
                to: 2048
            }

            STButton {
                enabled: !AppData.flux.ray_volume_flux_in_progress
                Layout.fillWidth: true
                Layout.columnSpan: 2

                text: "Start Raster"
                left_text_icon: "\uf0da"

                onClicked: {
                    AppData.flux.start_generate_volume_flux(resolution_spin.value)
                }
            }

            STPropertySeparator {
                title: "Isosurface"
            }

            STPropertyLabel {
                text: "Isovalue"
            }

            STDoubleSpinBox {
                id: iso_spin
                Layout.fillWidth: true

                value: 0.90
                from: 0.0
                stepSize: .01
                to: 1.0
            }

            STButton {
                enabled: !AppData.flux.ray_volume_flux_in_progress
                Layout.fillWidth: true
                Layout.columnSpan: 2

                text: "Generate Surface"
                left_text_icon: "\uf0da"

                onClicked: {
                    AppData.flux.start_generate_isosurface(iso_spin.value)
                }
            }

            STSwitch {
                text: "Visible"
                checked: AppData.flux.show_flux_volume
                onToggled: AppData.flux.show_flux_volume = checked
            }

        }
    }
}
