import QtCore
import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts
import QtQuick.Dialogs

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

        STFormPanel {
            collapsible: true
            title: "Current Flux Map"

            visible: AppData.flux.current_image.length || computed_map_combo.count > 0

            STFormRow {
                label: "Computed Map"
                visible: computed_map_combo.count > 0

                STComboBox {
                    id: computed_map_combo
                    Layout.fillWidth: true

                    model: AppData.flux.computed_maps_model
                    textRole: "name"
                    currentIndex: Math.max(
                                      0,
                                      AppData.flux.computed_maps_model.index_of(
                                          AppData.flux.current_entity))

                    onActivated: (index) => {
                        AppData.flux.select_entity(
                                    AppData.flux.computed_maps_model.entity_at(index))
                    }
                }
            }

            Image {
                id: flux_image_container
                Layout.fillWidth: true
                Layout.preferredHeight: width
                Layout.minimumHeight: 48

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

                    onClicked: save_image_dialog.open()

                    FileDialog {
                        id: save_image_dialog

                        fileMode: FileDialog.SaveFile
                        defaultSuffix: "png"
                        currentFolder: StandardPaths.standardLocations(
                                           StandardPaths.DocumentsLocation
                                           )[0]

                        onAccepted: {
                            AppData.flux.save_image(flux_image_container.source,
                                                    save_image_dialog.currentFile);
                        }
                    }
                }

                ShadowedRectangle {
                    width: 24
                    height: parent.height / 2
                    anchors.left: parent.left
                    anchors.bottom: parent.bottom

                    anchors.bottomMargin: 6
                    anchors.leftMargin: 6

                    opacity: legend_hover.containsMouse ? 1 : .5

                    Behavior on opacity {
                        NumberAnimation {
                            duration: 200
                        }
                    }

                    Image {
                        rotation: 90
                        source: "qrc:/assets/images/b_to_r_wide.png"
                        anchors.centerIn: parent
                        width: parent.height - 4
                        height: parent.width - 4
                    }

                    MouseArea {
                        id: legend_hover
                        hoverEnabled: true
                        anchors.fill: parent
                    }
                }
            }

            STComboBar {
                id: map_selector
                Layout.fillWidth: true
                iconModel: ["\uf00a", "\uf141"]
                model: ["Bins", "Points"]
            }

            STSwitch {
                text: "Show Whole Scene"
                checked: AppData.flux.show_other_geometry
                onToggled: AppData.flux.show_other_geometry = checked
            }

            STFormHelp {
                key: "analyze.flux.display"
            }

            Label {
                Layout.fillWidth: true
                text: "Flux map values are scaled by solar DNI when ray area metadata is available."
                color: App.theme.fontColor
                opacity: 0.75
                wrapMode: Label.WrapAtWordBoundaryOrAnywhere
            }

            STFormHelp {
                key: "analyze.flux.stats"
            }

            STFormRow {
                label: "Plotted Power (W)"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.plotted_power
                    font.bold: true
                }
            }

            STFormRow {
                label: "Peak Flux (W/m^2)"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.peak_flux
                    font.bold: true
                }
            }

            STFormRow {
                label: "Min Flux (W/m^2)"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.min_flux
                    font.bold: true
                }
            }

            STFormRow {
                label: "Average Flux (W/m^2)"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.average_flux
                    font.bold: true
                }
            }

            STFormRow {
                label: "Sigma Flux (W/m^2)"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.sigma_flux
                    font.bold: true
                }
            }

            STFormRow {
                label: "Uniformity"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.uniformity
                    font.bold: true
                }
            }

            STFormRow {
                label: "Peak Flux Uncert"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.peak_flux_uncertainty
                    font.bold: true
                }
            }

            STFormRow {
                label: "Average Flux Uncert"

                Label {
                    Layout.fillWidth: true
                    text: root.flux_module.current_flux_stats.average_flux_uncertainty
                    font.bold: true
                }
            }

            STFormRow {
                label: "Centroid"

                Label {
                    property vector3d cent: root.flux_module.current_flux_stats.centroid
                    Layout.fillWidth: true
                    text: root.formatCentroid(cent)
                    font.bold: true
                }
            }
        }

        STFormPanel {
            collapsible: true
            title: "Compute Flux Map"

            STFormRow {
                label: "Element"

                RowLayout {
                    Layout.fillWidth: true
                    spacing: 8

                    STButton {
                        Layout.fillWidth: true

                        text: AppData.flux.current_entity_name.length ?
                                  AppData.flux.current_entity_name : "Select Element"
                        left_text_icon: "\uf03a"

                        onClicked: entity_pop.open()

                        SelectItemPopup {
                            id: entity_pop
                            source_model: AppData.flux.entity_model

                            onSelectedEntity: (entity) => {
                                AppData.flux.select_entity(entity)
                            }
                        }
                    }

                    STIconButton {
                        icon: "\uf245"
                        toolTip: "Pick element from view"

                        onClicked: {
                            App.view.simulation_content_view = true
                            App.view.mouse_mode = ViewModule.SelectElement
                        }
                    }

                    STIconButton {
                        icon: "\uf140"
                        toolTip: "Orient Camera to Element"
                        enabled: AppData.flux.current_entity.is_valid()

                        onClicked: {
                            App.view.simulation_content_view = true
                            simulation_scene.orient_camera_to_database_position(
                                        AppData.flux.current_entity_position)
                        }
                    }
                }
            }

            STFormHelp {
                key: "analyze.flux.target"
            }

            STFormRow {
                label: "Solar DNI"

                STDoubleSpinBox {
                    Layout.fillWidth: true
                    from: 0.0
                    to: 2000.0
                    stepSize: 50.0
                    decimals: 1
                    value: AppData.flux.dni
                    suffix: " W/m^2"

                    onValueModified: AppData.flux.dni = value
                }
            }

            STFormHelp {
                key: "analyze.flux.dni"
            }

            STSwitch {
                text: "Show Triangle Grid"
                checked: AppData.flux.pending_flux_maps.show_mesh_grid
                onToggled: AppData.flux.pending_flux_maps.show_mesh_grid = checked
            }

            STFormRow {
                label: "Image Resolution (Pixels)"

                RowLayout {
                    Layout.fillWidth: true
                    spacing: 8

                    STSpinBox {
                        Layout.fillWidth: true
                        from: 64
                        to: 8192
                        stepSize: 64
                        value: AppData.flux.pending_flux_maps.image_resolution.width

                        onValueModified: {
                            AppData.flux.pending_flux_maps.image_resolution =
                                    Qt.size(value,
                                            AppData.flux.pending_flux_maps.image_resolution.height)
                        }
                    }

                    Label {
                        Layout.alignment: Qt.AlignVCenter
                        text: "\u00d7"
                        opacity: 0.7
                    }

                    STSpinBox {
                        Layout.fillWidth: true
                        from: 64
                        to: 8192
                        stepSize: 64
                        value: AppData.flux.pending_flux_maps.image_resolution.height

                        onValueModified: {
                            AppData.flux.pending_flux_maps.image_resolution =
                                    Qt.size(AppData.flux.pending_flux_maps.image_resolution.width,
                                            value)
                        }
                    }
                }
            }

            STFormRow {
                label: "Surface resolution"

                STSpinBox {
                    Layout.fillWidth: true
                    from: 1
                    to: 16
                    value: AppData.flux.pending_flux_maps.mesh_resolution_multiply

                    onValueModified: {
                        AppData.flux.pending_flux_maps.mesh_resolution_multiply = value
                    }
                }
            }

            STFormHelp {
                key: "analyze.flux.resolution"
            }

            ListView {
                Layout.fillWidth: true

                Layout.preferredHeight: 120

                model: root.flux_module.pending_flux_maps

                visible: count > 0

                delegate: Rectangle {
                    required property var entity;
                    required property int progress;

                    width: ListView.view.width
                    height: 24

                    color: Qt.alpha("black", .25)

                    radius: 10

                    border.color: Material.dividerColor
                    border.width: 1

                    ShadowedRectangle {
                        anchors.left: parent.left
                        anchors.top: parent.top
                        anchors.bottom: parent.bottom

                        width: progress / 100 * parent.width
                    }

                    Label {
                        anchors.fill: parent

                        text: "Map for " + parent.entity + ": " + parent.progress + "%"

                        verticalAlignment: Qt.AlignVCenter
                        horizontalAlignment: Qt.AlignHCenter
                    }

                }
            }

            STButton {
                Layout.fillWidth: true

                text: "Compute Map"
                left_text_icon: "\uf0da"

                onClicked: {
                    AppData.flux.start_generate()
                }
            }


        }

        STFormPanel {
            // Disable for now until this is a working feature
            visible: false

            collapsible: true
            title: "Normalized Ray Volume"

            STFormHelp {
                key: "analyze.flux.volume"
            }

            STFormRow {
                label: "Grid Resolution"

                STSpinBox {
                    id: resolution_spin
                    Layout.fillWidth: true

                    value: 512
                    from: 64
                    to: 2048
                }
            }

            STButton {
                enabled: !AppData.flux.ray_volume_flux_in_progress
                Layout.fillWidth: true

                text: "Start Raster"
                left_text_icon: "\uf0da"

                onClicked: {
                    AppData.flux.start_generate_volume_flux(resolution_spin.value)
                }
            }

            STPropertySeparator {
                title: "Isosurface"
            }

            STFormRow {
                label: "Isovalue"

                STDoubleSpinBox {
                    id: iso_spin
                    Layout.fillWidth: true

                    value: 0.90
                    from: 0.0
                    stepSize: .01
                    decimals: 2
                    to: 1.0
                }
            }

            STButton {
                enabled: !AppData.flux.ray_volume_flux_in_progress
                Layout.fillWidth: true

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
