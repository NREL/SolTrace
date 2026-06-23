import QtQuick
import QtQuick3D
import QtQuick3D.Helpers
import QtQuick3D.AssetUtils
import QtQuick.Controls.Material

import SolTrace

Node {
    id: world_node
    rotation: Quaternion.fromEulerAngles(-90, 0, 0)

    property real elevation: App.sun.type === SunModule.Directional ? sunDirectionRayGroup.z : pointSource.z

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

    Model {
        id: pointSource

        source: "#Sphere"

        x: App.sun.position.x
        y: App.sun.position.y
        z: App.sun.position.z
        scale: Qt.vector3d(App.view.sim.sun_viz_scale / 100, App.view.sim.sun_viz_scale / 100, App.view.sim.sun_viz_scale / 100)

        visible: App.sun.type === SunModule.PointSource && App.view.sim.sun_viz


        materials: [
            PrincipledMaterial {
                metalness: 1
                roughness: 0
                baseColor: App.view.sim.sun_color
            }
        ]
    }

    Node {
        id: sunDirectionRayGroup

        property vector3d sunDir: Qt.vector3d(App.sun.position.x,
                                              App.sun.position.y,
                                              App.sun.position.z)
        property int distance: 1000

        visible: App.sun.type === SunModule.Directional && App.view.sim.sun_viz
        scale: Qt.vector3d(App.view.sim.sun_viz_scale / 100, App.view.sim.sun_viz_scale / 100, App.view.sim.sun_viz_scale / 100)

        // Position the rays at the sun's location
        position: Qt.vector3d(sunDir.x * distance,
                              sunDir.y * distance,
                              sunDir.z * distance)

        // Orient so that local +Y points from sun toward origin.
        rotation: Quaternion.lookAt(
            position,                    // source: where we are
            Qt.vector3d(0, 0, 0),        // target: origin
            Qt.vector3d(0, 1, 0),        // forward = +Y (cylinder's long axis!)
            Qt.vector3d(0, 0, 1)         // up (any vector not parallel to forward)
        )

        // Ray arrows
        Repeater3D {
            model: parent.generatePositions()
            delegate: Node {
                required property var modelData
                position: modelData

                // Shaft
                Model {
                    source: "#Cylinder"
                    scale: Qt.vector3d(0.1, 5.0, 0.1)
                    materials: PrincipledMaterial {
                    metalness: 0
                    roughness: 1
                    baseColor: App.view.sim.sun_color
                    }
                }

                // Arrowhead
                Model {
                    source: "#Cone"
                    position: Qt.vector3d(0, 250, 0)
                    scale: Qt.vector3d(0.3, 0.5, 0.3)
                    materials: PrincipledMaterial {
                    metalness: 0
                    roughness: 1
                    baseColor: App.view.sim.sun_color
                    }
                }
            }
        }

        function generatePositions() {
            let rayPositions = []
            let rayGridSpan = 2
            let rayDistance = 100
            let rayGridOffset = rayDistance * (rayGridSpan - 1) / 2
            for (let i = 0; i < rayGridSpan; i++) {
                for (let j = 0; j < rayGridSpan; j++) {
                    let x = i * rayDistance - rayGridOffset
                    let z = j * rayDistance - rayGridOffset
                    rayPositions.push(Qt.vector3d(x, 0, z))
                }
            }
            return rayPositions
        }
    }
}