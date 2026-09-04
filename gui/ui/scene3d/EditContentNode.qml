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
        App.layout.world_geometry_model.set_surface_thickness(
                    App.view.sim.geometry_thickness)
        App.layout.world_geometry_model.set_subdivision_scale(
                    App.view.sim.geometry_subdivision_scale)
        App.layout.world_geometry_model.set_default_color(
                    App.view.sim.geometry_color)
    }

    Component.onCompleted: applyPerformanceSettings()

    Repeater3D {
        model: App.layout.world_geometry_model

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
            App.layout.world_geometry_model.set_default_color(App.view.sim.geometry_color)
        }
        function onGeometry_thickness_changed() {
            world_node.applyPerformanceSettings()
        }
        function onGeometry_subdivision_scale_changed() {
            world_node.applyPerformanceSettings()
        }
    }

    SunVisualizationNode {
        id: sunVisualization
        sourcePosition: Qt.vector3d(App.sun.position.x,
                                    App.sun.position.y,
                                    App.sun.position.z)
        isPointSource: App.sun.type === SunModule.PointSource

        Component.onCompleted: App.sun.update_position()
    }
}
