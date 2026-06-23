import QtQuick
import QtQuick3D
import QtQuick3D.Helpers

Node {
    id: root

    required property int activeAxis
    required property int gizmoMode

    property color xColor: "#e74c3c"
    property color xActiveColor: "#ff6b6b"
    property color yColor: "#3498db"
    property color yActiveColor: "#6bb5ff"
    property color zColor: "#2ecc71"
    property color zActiveColor: "#6bff9e"

    // Translation handles

    // X axis
    Model {
        source: "#Cylinder"
        objectName: "axis_0"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.04, 0.8, 0.04)
        position: Qt.vector3d(40, 0, 0)
        eulerRotation: Qt.vector3d(0, 0, -90)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 0 ? root.xActiveColor : root.xColor
            lighting: PrincipledMaterial.NoLighting
        }
    }
    Model {
        source: "#Cone"
        objectName: "axis_0"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.1, 0.2, 0.1)
        position: Qt.vector3d(70, 0, 0)
        eulerRotation: Qt.vector3d(0, 0, -90)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 0 ? root.xActiveColor : root.xColor
            lighting: PrincipledMaterial.NoLighting
        }
    }

    // Y axis
    Model {
        source: "#Cylinder"
        objectName: "axis_1"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.04, 0.8, 0.04)
        position: Qt.vector3d(0, 40, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 1 ? root.yActiveColor : root.yColor
            lighting: PrincipledMaterial.NoLighting
        }
    }
    Model {
        source: "#Cone"
        objectName: "axis_1"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.1, 0.2, 0.1)
        position: Qt.vector3d(0, 70, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 1 ? root.yActiveColor : root.yColor
            lighting: PrincipledMaterial.NoLighting
        }
    }

    // Z axis
    Model {
        source: "#Cylinder"
        objectName: "axis_2"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.04, 0.8, 0.04)
        position: Qt.vector3d(0, 0, 40)
        eulerRotation: Qt.vector3d(90, 0, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 2 ? root.zActiveColor : root.zColor
            lighting: PrincipledMaterial.NoLighting
        }
    }
    Model {
        source: "#Cone"
        objectName: "axis_2"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.1, 0.2, 0.1)
        position: Qt.vector3d(0, 0, 70)
        eulerRotation: Qt.vector3d(90, 0, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 2 ? root.zActiveColor : root.zColor
            lighting: PrincipledMaterial.NoLighting
        }
    }

    // XY plane handle
    Model {
        source: "#Rectangle"
        objectName: "plane_0"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.5, 0.5, 1)
        position: Qt.vector3d(30, 30, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 3 ? "#f5a623" : "#e67e22"
            lighting: PrincipledMaterial.NoLighting
            opacity: 0.5
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // XZ plane handle
    Model {
        source: "#Rectangle"
        objectName: "plane_1"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.5, 0.5, 1)
        position: Qt.vector3d(30, 0, 30)
        eulerRotation: Qt.vector3d(-90, 0, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 4 ? "#c39bd3" : "#9b59b6"
            lighting: PrincipledMaterial.NoLighting
            opacity: 0.5
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // YZ plane handle
    Model {
        source: "#Rectangle"
        objectName: "plane_2"
        pickable: true
        visible: gizmoMode === 0
        scale: Qt.vector3d(0.5, 0.5, 1)
        position: Qt.vector3d(0, 30, 30)
        eulerRotation: Qt.vector3d(0, 90, 0)
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 5 ? "#48d1b5" : "#1abc9c"
            lighting: PrincipledMaterial.NoLighting
            opacity: 0.5
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // X rotation ring
    Model {
        objectName: "rot_0"
        pickable: true
        visible: gizmoMode === 1
        eulerRotation: Qt.vector3d(0, 0, 90)
        geometry: TorusGeometry {
            radius: 65
            tubeRadius: 3
            rings: 48
            segments: 16
        }
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 0 ? root.xActiveColor : root.xColor
            lighting: PrincipledMaterial.NoLighting
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // Y rotation ring
    Model {
        objectName: "rot_1"
        pickable: true
        visible: gizmoMode === 1
        geometry: TorusGeometry {
            radius: 65
            tubeRadius: 3
            rings: 48
            segments: 16
        }
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 1 ? root.yActiveColor : root.yColor
            lighting: PrincipledMaterial.NoLighting
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // Z rotation ring
    Model {
        objectName: "rot_2"
        pickable: true
        visible: gizmoMode === 1
        eulerRotation: Qt.vector3d(90, 0, 0)
        geometry: TorusGeometry {
            radius: 65
            tubeRadius: 3
            rings: 48
            segments: 16
        }
        materials: PrincipledMaterial {
            baseColor: root.activeAxis === 2 ? root.zActiveColor : root.zColor
            lighting: PrincipledMaterial.NoLighting
            cullMode: PrincipledMaterial.NoCulling
        }
    }

    // Mode toggle sphere
    Model {
        source: "#Sphere"
        objectName: "mode_toggle"
        pickable: true
        scale: Qt.vector3d(0.12, 0.12, 0.12)
        materials: PrincipledMaterial {
            baseColor: root.gizmoMode === 0 ? "#cccccc" : "#ffcc00"
            lighting: PrincipledMaterial.NoLighting
        }
    }
}
