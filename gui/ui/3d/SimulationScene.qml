import QtQuick
import QtQuick3D
import QtQuick3D.Helpers
import QtQuick3D.AssetUtils

import SolTrace

Item {
    id: root

    property int activeAxis: -1
    property int gizmoMode: 0
    property bool showGizmo: App.view.mouse_mode === ViewModule.EditElement
                             && App.layout.instance_edit
    property bool isDragging: false
    property point lastMousePos: Qt.point(0, 0)
    property real initialAngle: 0.0
    property vector3d initialRotation: Qt.vector3d(0, 0, 0)

    readonly property var axisDirs: [
        Qt.vector3d(1, 0, 0),
        Qt.vector3d(0, 1, 0),
        Qt.vector3d(0, 0, 1)
    ]

    function align_to_axis(axis, invert) {
        controller.align_to_axis(axis, invert)
    }

    function scene_position_from_database_position(point) {
        return Qt.vector3d(point.x, point.z, -point.y)
    }

    function orient_camera_to_database_position(point) {
        controller.look_at(scene_position_from_database_position(point))
    }

    function reset_camera_view() {
        controller.reset_view()
    }

    View3D {
        id: view
        anchors.fill: parent
        camera: controller.active_camera

        environment: SceneEnvironment {
            antialiasingMode: SceneEnvironment.MSAA
            antialiasingQuality: SceneEnvironment.High

            temporalAAEnabled: true
            temporalAAStrength: 0.8

            aoStrength: 50
            aoDistance: 10
            aoSoftness: 75
            aoBias: 0.01
            aoSampleRate: 4

            probeExposure: 1.2

            tonemapMode: SceneEnvironment.TonemapModeAces

            backgroundMode: SceneEnvironment.SkyBox

            lightProbe: Texture {
                textureData: {
                    if (App.view.sim.sky === SimulationViewState.Day) return daySky
                    if (App.view.sim.sky === SimulationViewState.Blueprint) return blueprintSky
                    let elevation = App.view.simulation_content_view
                                  ? result_node.elevation
                                  : edit_node.elevation
                    if (elevation > 30) return daySky
                    if (elevation > 10) return lateAfternoonSky
                    if (elevation > -10) return sunsetSky
                    return nightSky
                }
                mappingMode: Texture.LightProbe
            }

            ProceduralSkyTextureData {
                id: daySky
                sunColor: Qt.rgba(0, 0, 0, 0)
                skyTopColor: Qt.rgba(0.2, 0.35, 0.6, 1.0)
                skyHorizonColor: Qt.rgba(0.55, 0.65, 0.75, 1.0)
                groundHorizonColor: Qt.rgba(0.55, 0.65, 0.75, 1.0)
                groundBottomColor: Qt.rgba(0.275, 0.325, 0.375, 1.0)
            }

            ProceduralSkyTextureData {
                id: lateAfternoonSky
                sunColor: Qt.rgba(0, 0, 0, 0)
                skyTopColor: Qt.rgba(0.3, 0.3, 0.5, 1.0)
                skyHorizonColor: Qt.rgba(0.75, 0.6, 0.5, 1.0)
                groundHorizonColor: Qt.rgba(0.45, 0.45, 0.55, 1.0)
                groundBottomColor: Qt.rgba(0.2, 0.2, 0.3, 1.0)
            }

            ProceduralSkyTextureData {
                id: sunsetSky
                sunColor: Qt.rgba(0, 0, 0, 0)
                skyTopColor: Qt.rgba(0.15, 0.15, 0.35, 1.0)
                skyHorizonColor: Qt.rgba(0.9, 0.5, 0.3, 1.0)
                groundHorizonColor: Qt.rgba(0.5, 0.35, 0.3, 1.0)
                groundBottomColor: Qt.rgba(0.15, 0.1, 0.15, 1.0)
            }

            ProceduralSkyTextureData {
                id: nightSky
                sunColor: Qt.rgba(0, 0, 0, 0)
                skyTopColor: Qt.rgba(0.02, 0.02, 0.08, 1.0)
                skyHorizonColor: Qt.rgba(0.05, 0.05, 0.15, 1.0)
                groundHorizonColor: Qt.rgba(0.05, 0.05, 0.1, 1.0)
                groundBottomColor: Qt.rgba(0.02, 0.02, 0.05, 1.0)
            }

            ProceduralSkyTextureData {
                id: blueprintSky
                sunColor: Qt.rgba(0, 0, 0, 0)
                skyTopColor: "#818182"
                skyHorizonColor: "#818182"
                groundHorizonColor: "#4d4d4d"
                groundBottomColor: "#4d4d4d"
            }

            InfiniteGrid {
                id: infiniteGrid
                visible: true
                gridInterval: 50
            }
        }

        PerspectiveCamera {
            id: mainPerspectiveCamera
            z: 100
            // TODO: Allow user customization
            clipNear: 1
        }

        OrthographicCamera {
            id: mainOrthoCamera
            z: 500
            clipNear: 0.01
        }

        DirectionalLight {
            eulerRotation.x: -45
            eulerRotation.y: 45
        }

        EditContentNode {
            id: edit_node
            visible: !App.view.simulation_content_view
        }

        SimulationResultNode {
            id: result_node
            visible: App.view.simulation_content_view
        }
    }

    GizmoOverlay {
        id: gizmoOverlay
        anchors.fill: parent
        sourceCamera: controller.active_camera
        gizmoPosition: App.layout.instance_edit ? App.layout.instance_edit.position : Qt.vector3d(0, 0, 0)
        activeAxis: root.activeAxis
        gizmoMode: root.gizmoMode
        enabled: root.showGizmo

        Binding {
            target: gizmoOverlay
            property: "enabled"
            value: root.showGizmo
        }
    }

    CameraController {
        id: controller

        enabled: !root.isDragging

        perspective_camera: mainPerspectiveCamera
        orthographic_camera: mainOrthoCamera
        default_perspective_position: Qt.vector3d(0, 0, 100)
        default_orthographic_position: Qt.vector3d(0, 0, 500)

        use_wasd: App.view.sim.camera === SimulationViewState.WASD
        use_orthographic: App.view.sim.perspective === SimulationViewState.Orthographic

        anchors.fill: parent

    }

    SimulationMouseArea {
        enabled: App.view.mouse_mode !== ViewModule.Camera
        view: view
        controller: controller
        gizmoOverlay: gizmoOverlay
        root: root
    }
}
