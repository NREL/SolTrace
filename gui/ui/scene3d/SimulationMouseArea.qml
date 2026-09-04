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
    property bool pendingScenePress: false
    property bool cameraDragActive: false
    property point pendingMousePos: Qt.point(0, 0)
    property bool pendingCameraPan: false
    property real cameraDragThreshold: 3.0

    SimulationPickController {
        id: pickController
        view: simMouseArea.view
        controller: simMouseArea.controller
        sceneState: simMouseArea.root
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
        var scenePos = pickController.toScene(worldPos)
        var sceneTip = pickController.toScene(Qt.vector3d(
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

    function traceGizmo(message) {
        console.debug("[SimulationMouseArea][gizmo]", message
                      + " isDragging=" + root.isDragging
                      + " activeAxis=" + root.activeAxis
                      + " mode=" + root.gizmoMode
                      + " mouseMode=" + App.view.mouse_mode)
    }

    function resetDeferredCameraDrag() {
        pendingScenePress = false
        cameraDragActive = false
        pendingCameraPan = false
    }

    function startDeferredCameraDrag(mouse) {
        cameraDragActive = true
        pendingScenePress = false
        controller.mousePressed(Qt.vector2d(pendingMousePos.x, pendingMousePos.y),
                                pendingCameraPan)
        controller.mouseMoved(Qt.vector2d(mouse.x, mouse.y))
        pickController.tracePick("deferred camera drag start x=" + mouse.x
                                 + " y=" + mouse.y
                                 + " pan=" + pendingCameraPan)
    }

    anchors.fill: parent
    acceptedButtons: App.view.mouse_mode === ViewModule.Camera
                     ? Qt.RightButton
                     : Qt.LeftButton | Qt.RightButton
    cursorShape: (App.view.mouse_mode === ViewModule.SelectElement
                  || App.view.mouse_mode === ViewModule.SelectMaterial
                  || App.view.mouse_mode === ViewModule.SelectGeometry
                  || App.view.mouse_mode === ViewModule.PickRay
                  || App.view.mouse_mode === ViewModule.SelectRayFilterElement)
                 ? Qt.CrossCursor
                 : Qt.ArrowCursor

    onPressed: (mouse) => {
        pickController.tracePick("pressed mode=" + App.view.mouse_mode
                                  + " button=" + mouse.button
                                  + " x=" + mouse.x + " y=" + mouse.y)

        if (mouse.button === Qt.RightButton) {
            pickController.selectEditableElementAt(mouse.x, mouse.y)
            return
        }

        if (App.view.mouse_mode === ViewModule.PickRay) {
            pickController.pickRay(mouse.x, mouse.y)
            pickController.returnToCameraModeIfOneShot()
            return
        }

        if (App.view.simulation_content_view
                && App.view.mouse_mode === ViewModule.SelectElement) {
            pickController.selectFluxElementFromResultView(mouse.x, mouse.y)
            return
        }

        if (App.view.simulation_content_view
                && App.view.mouse_mode === ViewModule.SelectRayFilterElement) {
            pickController.selectRayFilterElementFromResultView(mouse.x, mouse.y)
            return
        }

        // The analysis/simulation-result view has its own visual content. Other
        // mouse modes currently only support editing database geometry.
        if (App.view.simulation_content_view) {
            pickController.tracePick("ignored: simulation content view is active")
            pickController.returnToCameraModeIfOneShot()
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
                pickController.tracePick("gizmo hit objectName=" + name)

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
                    traceGizmo("begin axis drag objectName=" + name
                               + " x=" + mouse.x + " y=" + mouse.y)
                    return
                }

                // Translation within a plane. These are stored as activeAxis
                // 3..5 so onPositionChanged can distinguish plane movement
                // from single-axis movement.
                if (name.startsWith("plane_") && root.gizmoMode === 0) {
                    root.activeAxis = parseInt(name.split("_")[1]) + 3
                    root.isDragging = true
                    root.lastMousePos = Qt.point(mouse.x, mouse.y)
                    traceGizmo("begin plane drag objectName=" + name
                               + " x=" + mouse.x + " y=" + mouse.y)
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
                    traceGizmo("begin rotation drag objectName=" + name
                               + " x=" + mouse.x + " y=" + mouse.y
                               + " initialAngle=" + root.initialAngle
                               + " initialRotation=" + root.initialRotation)
                    return
                }
            }
        }

        if (App.view.mouse_mode === ViewModule.EditElement
                && mouse.button === Qt.LeftButton) {
            pendingScenePress = true
            pendingMousePos = Qt.point(mouse.x, mouse.y)
            pendingCameraPan =
                    Boolean(mouse.modifiers & Qt.ShiftModifier)
                    && !controller.use_wasd
            pickController.tracePick("deferred edit press x=" + mouse.x
                                      + " y=" + mouse.y
                                      + " pan=" + pendingCameraPan)
            return
        }

        pickController.handleScenePick(mouse.x, mouse.y, mouse.button)
    }

    onPositionChanged: (mouse) => {
        // Mouse motion only mutates geometry while a gizmo handle is actively
        // dragging. Ordinary mouse movement is left to CameraController.
        if (!root.isDragging || root.activeAxis < 0) {
            if (cameraDragActive) {
                controller.mouseMoved(Qt.vector2d(mouse.x, mouse.y))
                return
            }

            if (pendingScenePress) {
                var totalDx = mouse.x - pendingMousePos.x
                var totalDy = mouse.y - pendingMousePos.y
                var totalDistance =
                        Math.sqrt(totalDx * totalDx + totalDy * totalDy)
                if (totalDistance >= cameraDragThreshold) {
                    startDeferredCameraDrag(mouse)
                }
            }
            return
        }

        var ie = App.layout.instance_edit
        if (!ie) return

        var dx = mouse.x - root.lastMousePos.x
        var dy = mouse.y - root.lastMousePos.y

        traceGizmo("move x=" + mouse.x + " y=" + mouse.y
                   + " dx=" + dx + " dy=" + dy
                   + " buttons=" + mouse.buttons)

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
            traceGizmo("rotate currentAngle=" + currentAngle
                       + " deltaAngle=" + deltaAngle
                       + " result=" + ie.euler_angles_xyz)
            return
        }

        // Translation mode is incremental. Each move event consumes the delta
        // since the previous event and immediately updates lastMousePos.
        root.lastMousePos = Qt.point(mouse.x, mouse.y)

        // activeAxis 0..2: movement constrained to one world axis.
        if (root.activeAxis < 3) {
            var dir = axisDirs[root.activeAxis]
            var amount = projectMouseToAxis(dx, dy, dir)
            traceGizmo("translate axis=" + root.activeAxis
                       + " amount=" + amount
                       + " before=" + ie.position)
            ie.position = Qt.vector3d(
                ie.position.x + dir.x * amount,
                ie.position.y + dir.y * amount,
                ie.position.z + dir.z * amount
            )
            traceGizmo("translate after=" + ie.position)
        } else {
            // activeAxis 3..5: movement constrained to one of the principal
            // planes. Each plane is represented by a pair of world axes.
            var planeAxes = [[0, 1], [0, 2], [1, 2]]
            var axes = planeAxes[root.activeAxis - 3]
            var dir1 = axisDirs[axes[0]]
            var dir2 = axisDirs[axes[1]]
            var amount1 = projectMouseToAxis(dx, dy, dir1)
            var amount2 = projectMouseToAxis(dx, dy, dir2)
            traceGizmo("translate plane=" + root.activeAxis
                       + " amount1=" + amount1
                       + " amount2=" + amount2
                       + " before=" + ie.position)
            ie.position = Qt.vector3d(
                ie.position.x + dir1.x * amount1 + dir2.x * amount2,
                ie.position.y + dir1.y * amount1 + dir2.y * amount2,
                ie.position.z + dir1.z * amount1 + dir2.z * amount2
            )
            traceGizmo("translate after=" + ie.position)
        }
    }

    onReleased: (mouse) => {
        // Releasing the left button ends any gizmo drag and removes visual
        // active-axis highlighting.
        if (mouse.button === Qt.LeftButton) {
            if (cameraDragActive) {
                controller.mouseReleased(Qt.vector2d(mouse.x, mouse.y))
                pickController.tracePick("deferred camera drag release x="
                                          + mouse.x
                                          + " y=" + mouse.y)
                resetDeferredCameraDrag()
                return
            }

            if (pendingScenePress) {
                pickController.handleScenePick(mouse.x, mouse.y, mouse.button)
                resetDeferredCameraDrag()
                return
            }

            traceGizmo("release x=" + mouse.x + " y=" + mouse.y
                       + " buttons=" + mouse.buttons)
            root.isDragging = false
            root.activeAxis = -1
        }
    }

    onCanceled: {
        if (cameraDragActive) {
            controller.mouseReleased(Qt.vector2d(pendingMousePos.x,
                                                 pendingMousePos.y))
        }
        resetDeferredCameraDrag()
        traceGizmo("canceled")
        root.isDragging = false
        root.activeAxis = -1
    }
}
