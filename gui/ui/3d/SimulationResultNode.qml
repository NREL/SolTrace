import QtQuick
import QtQuick3D
import QtQuick3D.Helpers
import QtQuick3D.AssetUtils
import QtQuick.Controls.Material

import SolTrace


Node {
    id: world_node
    rotation: Quaternion.fromEulerAngles(-90, 0, 0)

    readonly property real elevation: sunVisualization.elevation
    readonly property bool blueprintMode: App.view.sim.sky === SimulationViewState.Blueprint

    function applyPerformanceSettings() {
        AppData.simulation.world_geometry_model.set_surface_thickness(
                    App.view.sim.geometry_thickness)
        AppData.simulation.world_geometry_model.set_subdivision_scale(
                    App.view.sim.geometry_subdivision_scale)
        AppData.simulation.world_geometry_model.set_default_color(
                    App.view.sim.geometry_color)
    }

    Component.onCompleted: applyPerformanceSettings()

    SunVisualizationNode {
        id: sunVisualization
        sourcePosition: App.simulation.result_sun_position
        isPointSource: App.simulation.result_sun_is_point_source
    }

    Repeater3D {
        visible: flux_repeater.count === 0 || AppData.flux.show_other_geometry
        model: AppData.simulation.world_geometry_model

        delegate: Model {
            id: geometry_model

            required property var group_instances
            required property var group_geometry
            property bool is_focused: false

            instancing: group_instances
            geometry: group_geometry
            pickable: true

            materials: [
                PrincipledMaterial {
                    metalness: world_node.blueprintMode ? 0 : 1
                    roughness: world_node.blueprintMode ? 1 : 0
                    baseColor: "white"

                    lighting: world_node.blueprintMode ? PrincipledMaterial.NoLighting : PrincipledMaterial.FragmentLighting
                }
            ]
        }
    }

    Connections {
        target: App.view.sim
        function onGeometry_color_changed() {
            AppData.simulation.world_geometry_model.set_default_color(App.view.sim.geometry_color)
        }
        function onGeometry_thickness_changed() {
            world_node.applyPerformanceSettings()
        }
        function onGeometry_subdivision_scale_changed() {
            world_node.applyPerformanceSettings()
        }
    }

    Model {
        visible: AppData.view.show_intersections
        geometry: AppData.intersections.ray_geometry

        materials : [
            PrincipledMaterial {
                metalness: 0
                roughness: 1
                lighting: PrincipledMaterial.NoLighting
                baseColor: "white"
                baseColorMap: Texture {
                    source: "qrc:/assets/images/b_to_r_wide.png"
                }
            }
        ]
    }

    Repeater3D {
        id: flux_repeater
        model: AppData.flux.flux_map_world_model

        delegate: Model {
            required property var flux_position
            required property var flux_rotation
            required property var flux_geometry
            required property var flux_texture_data
            required property string flux_image_path

            geometry: flux_geometry
            //source: "#Cube"

            position: flux_position
            rotation: flux_rotation

            Component.onCompleted: {
                AppData.flux.current_image = flux_image_path
            }

            materials: [
                PrincipledMaterial {
                    metalness: 1
                    roughness: 0
                    baseColor: "white"

                    lighting: PrincipledMaterial.NoLighting

                    cullMode: PrincipledMaterial.NoCulling
                    baseColorMap: Texture {
                        textureData: flux_texture_data
                    }
                }
            ]
        }
    }

    Model {
        id: iso_vol_mesh

        visible: AppData.flux.show_flux_volume

        geometry: AppData.flux.ray_iso_volume

        materials: [
            PrincipledMaterial {
                metalness: 0
                roughness: .5
                baseColor: "white"
                cullMode: PrincipledMaterial.NoCulling
            }
        ]
    }

}
