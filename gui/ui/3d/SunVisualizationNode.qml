import QtQuick
import QtQuick3D

import SolTrace

Node {
    id: root

    property vector3d sourcePosition: Qt.vector3d(0, 0, 1)
    property bool isPointSource: false

    readonly property real sourceLength: Math.sqrt(
                                             sourcePosition.x * sourcePosition.x
                                             + sourcePosition.y * sourcePosition.y
                                             + sourcePosition.z * sourcePosition.z)
    readonly property vector3d sourceDirection: sourceLength > 1.0e-12
                                                  ? Qt.vector3d(sourcePosition.x / sourceLength,
                                                               sourcePosition.y / sourceLength,
                                                               sourcePosition.z / sourceLength)
                                                  : Qt.vector3d(0, 0, 1)

    readonly property real elevation: isPointSource
                                      ? pointSourceModel.z
                                      : directionalGroup.z

    Model {
        id: pointSourceModel
        source: "#Sphere"
        position: root.sourcePosition
        scale: Qt.vector3d(App.view.sim.sun_viz_scale / 100,
                           App.view.sim.sun_viz_scale / 100,
                           App.view.sim.sun_viz_scale / 100)
        visible: root.isPointSource && App.view.sim.sun_viz

        materials: PrincipledMaterial {
            metalness: 0
            roughness: 1
            baseColor: App.view.sim.sun_color
        }
    }

    Node {
        id: directionalGroup

        property int distance: 1000

        visible: !root.isPointSource && App.view.sim.sun_viz
        scale: Qt.vector3d(App.view.sim.sun_viz_scale / 100,
                           App.view.sim.sun_viz_scale / 100,
                           App.view.sim.sun_viz_scale / 100)
        position: Qt.vector3d(root.sourceDirection.x * distance,
                              root.sourceDirection.y * distance,
                              root.sourceDirection.z * distance)

        // The arrow geometry points along local +Y.
        rotation: Quaternion.lookAt(position,
                                    Qt.vector3d(0, 0, 0),
                                    Qt.vector3d(0, 1, 0),
                                    Qt.vector3d(0, 0, 1))

        Repeater3D {
            model: directionalGroup.generatePositions()

            delegate: Node {
                required property var modelData
                position: modelData

                Model {
                    source: "#Cylinder"
                    scale: Qt.vector3d(0.1, 5.0, 0.1)
                    materials: PrincipledMaterial {
                        metalness: 0
                        roughness: 1
                        baseColor: App.view.sim.sun_color
                    }
                }

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
