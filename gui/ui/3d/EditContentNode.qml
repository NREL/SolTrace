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
                    metalness: App.view.sim.blueprint_mode ? 0 : 1
                    roughness: App.view.sim.blueprint_mode ? 1 : 0
                    baseColor: App.view.sim.geometry_color

                    lighting: App.view.sim.blueprint_mode ? PrincipledMaterial.NoLighting : PrincipledMaterial.FragmentLighting
                }
            ]
        }
    }

    Connections {
        target: App.view.sim
        function onGeometry_color_changed() {
            App.layout.world_geometry_model.set_all_color(App.view.sim.geometry_color)
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
