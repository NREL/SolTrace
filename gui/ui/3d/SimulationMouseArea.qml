import QtQuick

// Owns mouse interaction for the editable 3D scene.
//
// This component deliberately sits outside the View3D scene graph. It receives
// 2D mouse events, asks View3D/GizmoOverlay what 3D object is under the cursor,
// then updates application state:
// - clicking a scene instance selects it for layout editing
// - clicking a scene instance in material/geometry mode selects that group
// - clicking in ray-pick mode sends a camera ray to the intersections backend
// - dragging a transform gizmo axis/plane changes the edited instance
//
// Coordinate note:
// The simulation/database uses x/y/z with z as vertical. The Qt Quick 3D scene
// is rotated in EditContentNode, so helper functions convert selected object
// positions into scene-space before asking Qt Quick 3D to project them.
MouseArea {
    id: simMouseArea

    // View3D used for scene picking and viewport dimensions.
    required property var view

    // CameraController exposes the currently active camera for projection.
    required property var controller

    // Separate overlay View3D containing only transform gizmo geometry.
    // Gizmo picking is done before scene picking so handles win over meshes.
    required property var gizmoOverlay

    // SimulationScene root object. Shared interaction state such as activeAxis,
    // gizmoMode, and isDragging lives there because the visible gizmo also binds
    // to it.
    required property var root

    // SolTrace/database-space unit axes. activeAxis 0..2 refers to these axes.
    readonly property var axisDirs: [
        Qt.vector3d(1, 0, 0),
        Qt.vector3d(0, 1, 0),
        Qt.vector3d(0, 0, 1)
    ]

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

    // Convert a 2D mouse delta into movement along a 3D world axis.
    //
    // Procedure:
    // 1. Project the selected object's position to screen space.
    // 2. Project a second point 100 units along the requested axis.
    // 3. The screen-space vector between those points tells us how that axis
    //    appears to the user from the current camera.
    // 4. Dot the mouse delta with that screen-space axis vector.
    //
    // The result is intentionally approximate. It makes drag direction feel
    // aligned with the handle on screen; it is not ray/plane intersection math.
    function projectMouseToAxis(dx, dy, axisDir) {
        var ie = App.layout.instance_edit
        if (!ie) return 0
        var worldPos = ie.position
        var cam = controller.active_camera
        var scenePos = toScene(worldPos)
        var sceneTip = toScene(Qt.vector3d(
            worldPos.x + axisDir.x * 100,
            worldPos.y + axisDir.y * 100,
            worldPos.z + axisDir.z * 100
        ))
        var screenOrigin = cam.mapToViewport(scenePos)
        var screenTip = cam.mapToViewport(sceneTip)
        var screenDirX = (screenTip.x - screenOrigin.x) * view.width
        var screenDirY = (screenTip.y - screenOrigin.y) * view.height
        var len = Math.sqrt(screenDirX * screenDirX + screenDirY * screenDirY)
        if (len < 0.001) return 0
        screenDirX /= len
        screenDirY /= len
        return (dx * screenDirX + dy * screenDirY) * 0.5
    }

    // Return the cursor angle around the selected object in screen space.
    // Rotation handles use this angle difference rather than x/y deltas so a
    // circular mouse motion maps naturally to a rotation amount.
    function screenAngleToObject(mx, my) {
        var ie = App.layout.instance_edit
        if (!ie) return 0
        var cam = controller.active_camera
        var sp = cam.mapToViewport(ie.position)
        var cx = sp.x * view.width
        var cy = sp.y * view.height
        return Math.atan2(my - cy, mx - cx) * (180.0 / Math.PI)
    }

    function returnToCameraModeIfOneShot() {
        if (App.view.mouse_mode === ViewModule.SelectElement
                || App.view.mouse_mode === ViewModule.SelectMaterial
                || App.view.mouse_mode === ViewModule.SelectGeometry
                || App.view.mouse_mode === ViewModule.PickRay) {
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

    anchors.fill: parent
    acceptedButtons: Qt.LeftButton
    cursorShape: (App.view.mouse_mode === ViewModule.SelectElement
                  || App.view.mouse_mode === ViewModule.SelectMaterial
                  || App.view.mouse_mode === ViewModule.SelectGeometry
                  || App.view.mouse_mode === ViewModule.PickRay)
                 ? Qt.CrossCursor
                 : Qt.ArrowCursor

    onPressed: (mouse) => {
        tracePick("pressed mode=" + App.view.mouse_mode
                  + " x=" + mouse.x + " y=" + mouse.y)

        if (App.view.mouse_mode === ViewModule.PickRay) {
            pickRay(mouse.x, mouse.y)
            returnToCameraModeIfOneShot()
            return
        }

        if (App.view.simulation_content_view
                && App.view.mouse_mode === ViewModule.SelectElement) {
            selectFluxElementFromResultView(mouse.x, mouse.y)
            return
        }

        // The analysis/simulation-result view has its own visual content. Other
        // mouse modes currently only support editing database geometry.
        if (App.view.simulation_content_view) {
            tracePick("ignored: simulation content view is active")
            returnToCameraModeIfOneShot()
            return
        }

        // First priority: transform gizmo handles. When a handle is hit, this
        // starts a drag operation and records enough state to interpret later
        // onPositionChanged events.
        if (App.view.mouse_mode === ViewModule.EditElement
                && root.showGizmo
                && mouse.button === Qt.LeftButton) {
            var gizmoResult = gizmoOverlay.pick(mouse.x, mouse.y)
            if (gizmoResult.objectHit) {
                var name = gizmoResult.objectHit.objectName
                tracePick("gizmo hit objectName=" + name)

                // Center handle switches between translate and rotate modes.
                if (name === "mode_toggle") {
                    root.gizmoMode = (root.gizmoMode === 0) ? 1 : 0
                    root.activeAxis = -1
                    return
                }

                // Translation along one axis. TransformGizmo names these
                // objects axis_0, axis_1, axis_2.
                if (name.startsWith("axis_") && root.gizmoMode === 0) {
                    root.activeAxis = parseInt(name.split("_")[1])
                    root.isDragging = true
                    root.lastMousePos = Qt.point(mouse.x, mouse.y)
                    return
                }

                // Translation within a plane. These are stored as activeAxis
                // 3..5 so onPositionChanged can distinguish plane movement
                // from single-axis movement.
                if (name.startsWith("plane_") && root.gizmoMode === 0) {
                    root.activeAxis = parseInt(name.split("_")[1]) + 3
                    root.isDragging = true
                    root.lastMousePos = Qt.point(mouse.x, mouse.y)
                    return
                }

                // Rotation around one axis. We keep the initial screen angle
                // and selected object's Euler rotation so dragging can apply a
                // delta relative to the starting pose.
                if (name.startsWith("rot_") && root.gizmoMode === 1) {
                    root.activeAxis = parseInt(name.split("_")[1])
                    root.isDragging = true
                    root.lastMousePos = Qt.point(mouse.x, mouse.y)
                    root.initialAngle = screenAngleToObject(mouse.x, mouse.y)
                    var ie = App.layout.instance_edit
                    if (ie) {
                        root.initialRotation = ie.euler_angles_xyz
                    }
                    return
                }
            }
        }

        // Second priority: editable scene geometry. Pick returns the Model that
        // was hit plus an instanceIndex. For instanced geometry the model is a
        // geometry group and instanceIndex selects the actual element instance.
        const result = view.pick(mouse.x, mouse.y)
        var object = result.objectHit
        if (!object) {
            tracePick("scene pick miss")
            root.activeAxis = -1
            returnToCameraModeIfOneShot()
            return
        }

        tracePick("scene pick hit object=" + object
                  + " hasInstancing=" + Boolean(object.instancing)
                  + " instanceIndex=" + result.instanceIndex)

        if (!object.instancing && mouse.button === Qt.LeftButton) {
            tracePick("scene pick hit non-instanced object")
            returnToCameraModeIfOneShot()
        } else if (object.instancing && mouse.button === Qt.LeftButton) {
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
                // Selecting a concrete geometry instance makes it the layout
                // edit target and opens the left panel directly to Layout editing.
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

    onPositionChanged: (mouse) => {
        // Mouse motion only mutates geometry while a gizmo handle is actively
        // dragging. Ordinary mouse movement is left to CameraController.
        if (!root.isDragging || root.activeAxis < 0) return

        var ie = App.layout.instance_edit
        if (!ie) return

        var dx = mouse.x - root.lastMousePos.x
        var dy = mouse.y - root.lastMousePos.y

        // Rotation mode uses total angular difference from the drag start. This
        // avoids accumulating rounding error and lets the object follow circular
        // cursor motion around its projected center.
        if (root.gizmoMode === 1 && root.activeAxis >= 0 && root.activeAxis < 3) {
            var currentAngle = screenAngleToObject(mouse.x, mouse.y)
            var deltaAngle = currentAngle - root.initialAngle

            // Keep the delta continuous when the cursor crosses +/-180 degrees.
            while (deltaAngle > 180) deltaAngle -= 360
            while (deltaAngle < -180) deltaAngle += 360

            // Pick a sign based on which side of the object the camera is on.
            // This keeps clockwise/counter-clockwise dragging feeling stable
            // when looking from opposite sides of an axis.
            var camPos = controller.active_camera.position
            var objPos = ie.position
            var toCamera = Qt.vector3d(
                camPos.x - objPos.x,
                camPos.y - objPos.y,
                camPos.z - objPos.z
            )

            var sign = 1.0
            if (root.activeAxis === 0) sign = toCamera.x >= 0 ? -1.0 : 1.0
            else if (root.activeAxis === 1) sign = toCamera.y >= 0 ? -1.0 : 1.0
            else if (root.activeAxis === 2) sign = toCamera.z >= 0 ? -1.0 : 1.0

            var rx = root.initialRotation.x
            var ry = root.initialRotation.y
            var rz = root.initialRotation.z

            if (root.activeAxis === 0) rx += deltaAngle * sign
            else if (root.activeAxis === 1) ry += deltaAngle * sign
            else if (root.activeAxis === 2) rz += deltaAngle * sign

            ie.euler_angles_xyz = Qt.vector3d(rx, ry, rz)
            return
        }

        // Translation mode is incremental. Each move event consumes the delta
        // since the previous event and immediately updates lastMousePos.
        root.lastMousePos = Qt.point(mouse.x, mouse.y)

        // activeAxis 0..2: movement constrained to one world axis.
        if (root.activeAxis < 3) {
            var dir = axisDirs[root.activeAxis]
            var amount = projectMouseToAxis(dx, dy, dir)
            ie.position = Qt.vector3d(
                ie.position.x + dir.x * amount,
                ie.position.y + dir.y * amount,
                ie.position.z + dir.z * amount
            )
        } else {
            // activeAxis 3..5: movement constrained to one of the principal
            // planes. Each plane is represented by a pair of world axes.
            var planeAxes = [[0, 1], [0, 2], [1, 2]]
            var axes = planeAxes[root.activeAxis - 3]
            var dir1 = axisDirs[axes[0]]
            var dir2 = axisDirs[axes[1]]
            var amount1 = projectMouseToAxis(dx, dy, dir1)
            var amount2 = projectMouseToAxis(dx, dy, dir2)
            ie.position = Qt.vector3d(
                ie.position.x + dir1.x * amount1 + dir2.x * amount2,
                ie.position.y + dir1.y * amount1 + dir2.y * amount2,
                ie.position.z + dir1.z * amount1 + dir2.z * amount2
            )
        }
    }

    onReleased: (mouse) => {
        // Releasing the left button ends any gizmo drag and removes visual
        // active-axis highlighting.
        if (mouse.button === Qt.LeftButton) {
            root.isDragging = false
            root.activeAxis = -1
        }
    }
}
