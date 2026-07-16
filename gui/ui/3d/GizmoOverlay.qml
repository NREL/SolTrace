import QtQuick
import QtQuick3D
import SolTrace

View3D {
    id: root

    required property var sourceCamera
    required property vector3d gizmoPosition
    required property int activeAxis
    required property int gizmoMode
    required property bool enabled

    environment: SceneEnvironment {
        backgroundMode: SceneEnvironment.Transparent
    }

    PerspectiveCamera {
        id: gizmoCamera
        position: root.sourceCamera ? root.sourceCamera.position : Qt.vector3d(0, 0, 0)
        eulerRotation: root.sourceCamera ? root.sourceCamera.eulerRotation : Qt.vector3d(0, 0, 0)
        fieldOfView: root.sourceCamera ? root.sourceCamera.fieldOfView : 60
        clipNear: root.sourceCamera ? root.sourceCamera.clipNear : 1
        clipFar: root.sourceCamera ? root.sourceCamera.clipFar : 10000
    }

    TransformGizmo {
        position: Qt.vector3d(root.gizmoPosition.x, root.gizmoPosition.z, -root.gizmoPosition.y)
        eulerRotation: Qt.vector3d(-90, 0, 0)
        visible: root.enabled
        activeAxis: root.activeAxis
        gizmoMode: root.gizmoMode
        scale: Qt.vector3d(0.5, 0.5, 0.5)
    }
}
