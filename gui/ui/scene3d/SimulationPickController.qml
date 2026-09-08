import QtQuick

import SolTrace

QtObject {
    id: root

    required property var view
    required property var controller
    required property var sceneState

    // Convert a point from database/world coordinates into the coordinate space
    // used by the Qt Quick 3D scene. This mirrors EditContentNode's -90 degree
    // rotation and is only used for camera projection math.
    function toScene(p) {
        return Qt.vector3d(p.x, p.z, -p.y)
    }

    // Inverse of toScene(). Use this when passing a camera ray to C++ result
    // data, which is stored in SolTrace/database coordinates.
    function fromScene(p) {
        return Qt.vector3d(p.x, -p.z, p.y)
    }

    function returnToCameraModeIfOneShot() {
        if (App.view.mouse_mode === ViewModule.SelectElement
                || App.view.mouse_mode === ViewModule.SelectMaterial
                || App.view.mouse_mode === ViewModule.SelectGeometry
                || App.view.mouse_mode === ViewModule.PickRay
                || App.view.mouse_mode === ViewModule.SelectRayFilterElement) {
            App.view.mouse_mode = ViewModule.Camera
        }
    }

    function tracePick(message) {
        console.log("[SimulationMouseArea]", message)
    }

    function entityString(entity) {
        if (!entity) return "<null>"
        if (entity.debug_string) return entity.debug_string()
        if (entity.value !== undefined) return "entity(" + entity.value + ")"
        return String(entity)
    }

    function clickRay(mx, my) {
        var cam = controller.active_camera
        var originScene = Qt.vector3d(cam.scenePosition.x,
                                      cam.scenePosition.y,
                                      cam.scenePosition.z)
        var forward = Qt.vector3d(cam.forward.x,
                                  cam.forward.y,
                                  cam.forward.z).normalized()
        var right = Qt.vector3d(cam.right.x,
                                cam.right.y,
                                cam.right.z).normalized()
        var up = Qt.vector3d(cam.up.x,
                             cam.up.y,
                             cam.up.z).normalized()

        var nx = (mx / Math.max(1, view.width) - 0.5) * 2.0
        var ny = (0.5 - my / Math.max(1, view.height)) * 2.0

        var fovDegrees = cam.fieldOfView !== undefined ? cam.fieldOfView : 45.0
        var fov = fovDegrees * Math.PI / 180.0
        var halfHeight = Math.tan(fov * 0.5)
        var halfWidth = halfHeight * Math.max(1, view.width) / Math.max(1, view.height)

        var directionScene = forward
            .plus(right.times(nx * halfWidth))
            .plus(up.times(ny * halfHeight))
            .normalized()

        return {
            position: fromScene(originScene),
            direction: fromScene(directionScene).normalized()
        }
    }

    function pickRay(mx, my) {
        var ray = clickRay(mx, my)
        tracePick("pick ray position=" + ray.position
                  + " direction=" + ray.direction)
        AppData.intersections.ray_geometry.pick_ray(ray.position,
                                                    ray.direction,
                                                    0.01)
    }

    function openLayoutEditorFor(entity) {
        tracePick("select element -> " + entityString(entity))
        App.view.workflow_phase = ViewModule.Configure
        App.view.left_panel.visible = true
        App.view.configure_section = 3
        App.view.editing_layout = true
        App.layout.edited_element = entity
    }

    function openMaterialEditorFor(entity) {
        tracePick("select material -> " + entityString(entity))
        App.view.workflow_phase = ViewModule.Configure
        App.view.left_panel.visible = true
        App.view.configure_section = 1
        App.view.editing_material = true
        App.materials.current_material = entity
    }

    function openGeometryEditorFor(entity) {
        tracePick("select geometry -> " + entityString(entity))
        App.view.workflow_phase = ViewModule.Configure
        App.view.left_panel.visible = true
        App.view.configure_section = 2
        App.view.editing_geometry = true
        App.materials.current_geometry = entity
    }

    function selectFluxElementFromResultView(mx, my) {
        const result = view.pick(mx, my)
        var object = result.objectHit
        if (!object || !object.instancing) {
            tracePick("flux element pick miss")
            returnToCameraModeIfOneShot()
            return
        }

        const index = result.instanceIndex
        if (index < 0) {
            tracePick("flux element pick failed: invalid instanceIndex=" + index)
            returnToCameraModeIfOneShot()
            return
        }

        var elementEntity = object.instancing.at(index)
        tracePick("flux element pick -> " + entityString(elementEntity))
        AppData.flux.select_entity(elementEntity)
        App.view.workflow_phase = ViewModule.Analyze
        App.view.left_panel.visible = true
        returnToCameraModeIfOneShot()
    }

    function selectRayFilterElementFromResultView(mx, my) {
        const result = view.pick(mx, my)
        var object = result.objectHit
        if (!object || !object.instancing) {
            tracePick("ray filter element pick miss")
            returnToCameraModeIfOneShot()
            return
        }

        const index = result.instanceIndex
        if (index < 0) {
            tracePick("ray filter element pick failed: invalid instanceIndex=" + index)
            returnToCameraModeIfOneShot()
            return
        }

        var elementEntity = object.instancing.at(index)
        tracePick("ray filter element pick -> " + entityString(elementEntity))
        AppData.intersections.ray_geometry.select_entity_filter(elementEntity)
        App.view.workflow_phase = ViewModule.Analyze
        App.view.left_panel.visible = true
        returnToCameraModeIfOneShot()
    }

    function handleScenePick(mx, my, button) {
        const result = view.pick(mx, my)
        var object = result.objectHit
        if (!object) {
            tracePick("scene pick miss")
            sceneState.activeAxis = -1
            returnToCameraModeIfOneShot()
            return
        }

        tracePick("scene pick hit object=" + object
                  + " hasInstancing=" + Boolean(object.instancing)
                  + " instanceIndex=" + result.instanceIndex)

        if (!object.instancing && button === Qt.LeftButton) {
            tracePick("scene pick hit non-instanced object")
            returnToCameraModeIfOneShot()
        } else if (object.instancing && button === Qt.LeftButton) {
            const index = result.instanceIndex
            if (index < 0) {
                tracePick("scene pick failed: invalid instanceIndex=" + index)
                returnToCameraModeIfOneShot()
                return
            }

            var elementEntity = object.instancing.at(index)
            var materialEntity = object.instancing.material_of(index)
            var geometryEntity = object.instancing.geometry_of(index)

            tracePick("instance index=" + index
                      + " element=" + entityString(elementEntity)
                      + " material=" + entityString(materialEntity)
                      + " geometry=" + entityString(geometryEntity))

            if (App.view.mouse_mode === ViewModule.SelectElement
                    || App.view.mouse_mode === ViewModule.EditElement) {
                openLayoutEditorFor(elementEntity)
            } else if (App.view.mouse_mode === ViewModule.SelectMaterial) {
                openMaterialEditorFor(materialEntity)
            } else if (App.view.mouse_mode === ViewModule.SelectGeometry) {
                openGeometryEditorFor(geometryEntity)
            } else {
                tracePick("scene pick ignored for mode=" + App.view.mouse_mode)
            }

            returnToCameraModeIfOneShot()
        }
    }

    function selectEditableElementAt(mx, my) {
        if (App.view.simulation_content_view) {
            tracePick("right-click select ignored: simulation content view is active")
            return
        }

        const result = view.pick(mx, my)
        var object = result.objectHit
        if (!object || !object.instancing) {
            tracePick("right-click select miss")
            return
        }

        const index = result.instanceIndex
        if (index < 0) {
            tracePick("right-click select failed: invalid instanceIndex=" + index)
            return
        }

        var elementEntity = object.instancing.at(index)
        tracePick("right-click select -> " + entityString(elementEntity))
        openLayoutEditorFor(elementEntity)
    }
}
