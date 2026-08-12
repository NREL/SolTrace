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

    function empty_bounds() {
        return {
            valid: false,
            min: Qt.vector3d(0, 0, 0),
            max: Qt.vector3d(0, 0, 0)
        }
    }

    function include_scene_point(bounds, point) {
        if (!bounds.valid) {
            bounds.min = Qt.vector3d(point.x, point.y, point.z)
            bounds.max = Qt.vector3d(point.x, point.y, point.z)
            bounds.valid = true
            return
        }

        bounds.min = Qt.vector3d(Math.min(bounds.min.x, point.x),
                                 Math.min(bounds.min.y, point.y),
                                 Math.min(bounds.min.z, point.z))
        bounds.max = Qt.vector3d(Math.max(bounds.max.x, point.x),
                                 Math.max(bounds.max.y, point.y),
                                 Math.max(bounds.max.z, point.z))
    }

    function include_database_bounds(bounds, source) {
        if (!source || !source.valid)
            return

        var minPoint = source.min
        var maxPoint = source.max

        for (var ix = 0; ix < 2; ++ix) {
            for (var iy = 0; iy < 2; ++iy) {
                for (var iz = 0; iz < 2; ++iz) {
                    include_scene_point(bounds, scene_position_from_database_position(
                                            Qt.vector3d(ix ? maxPoint.x : minPoint.x,
                                                        iy ? maxPoint.y : minPoint.y,
                                                        iz ? maxPoint.z : minPoint.z)))
                }
            }
        }
    }

    function fit_all_in_view() {
        var bounds = empty_bounds()

        if (App.view.simulation_content_view) {
            var fluxBounds = AppData.flux.flux_map_world_model.content_bounds()
            var normalGeometryVisible = AppData.flux.show_other_geometry
                    || !fluxBounds.valid

            if (normalGeometryVisible) {
                include_database_bounds(
                            bounds,
                            AppData.simulation.world_geometry_model.content_bounds(false))
            }

            include_database_bounds(bounds, fluxBounds)

            if (AppData.view.show_intersections) {
                include_database_bounds(bounds,
                                        AppData.simulation.current_result_bounds())
            }
        } else {
            include_database_bounds(
                        bounds,
                        AppData.layout.world_geometry_model.content_bounds(false))
        }

        if (!bounds.valid) {
            reset_camera_view()
            return
        }

        controller.fit_bounds(bounds.min, bounds.max)
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

            clearColor: "#0041BA"

            tonemapMode: SceneEnvironment.TonemapModeAces

            backgroundMode: {
                if (App.view.sim.sky === SimulationViewState.Blueprint) {
                    return SceneEnvironment.Color
                }

                return SceneEnvironment.SkyBox
            }

            lightProbe: App.view.sim.sky === SimulationViewState.Realistic
                        ? hdriSky
                        : proceduralSky

            Texture {
                id: proceduralSky
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

            Texture {
                id: hdriSky
                source: "qrc:/assets/skyboxes/" + ["clear_puresky", "partly_cloudy", "sunset_puresky"][App.view.sim.realistic_sky] + "_1k_sand.hdr"
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
                visible: App.view.sim.show_grid
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
            z: 1000
            clipNear: 0.01
            horizontalMagnification: 100
            verticalMagnification: 100
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
        input_enabled: !App.view.full_panel.visible

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
