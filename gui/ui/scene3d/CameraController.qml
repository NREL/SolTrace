import QtQuick
import QtQuick3D
import SolTrace

// CameraController owns the user-facing camera navigation modes used by the
// 3D views. It intentionally contains both navigation models in one place so
// the scene can swap between them without changing the cameras themselves:
//
//   - Orbit mode: drag to orbit around a target point, Shift-drag to pan, and
//     wheel to dolly toward/away from the target.
//   - WASD mode: drag to look around, keyboard to fly, and wheel to adjust
//     movement speed.
//
// The controller does not render anything. It is an invisible input surface
// that mutates the camera nodes passed in by the owning scene.
Item {
    id: root

    enum Axis {
        X,
        Y,
        Z
    }

    // Perspective camera is required because every controller instance needs at
    // least one concrete camera to drive.
    required property PerspectiveCamera perspective_camera

    // Orthographic camera is optional. When present, use_orthographic selects
    // it as active_camera. Orthographic mode is kept orbit-only because WASD
    // translation semantics are confusing without perspective depth.
    property OrthographicCamera orthographic_camera

    readonly property Camera active_camera : use_orthographic ? orthographic_camera : perspective_camera

    // Orbit mode uses rotation_target as the base point of interest. If no
    // target is supplied, the world origin is used.
    property Node rotation_target

    // Panning does not move rotation_target itself. Instead, it accumulates an
    // offset so pan can be layered on top of whatever scene node is currently
    // acting as the target. This matters for previews where rotation_target is
    // a stable helper Node owned elsewhere.
    property vector3d rotation_target_offset: Qt.vector3d(0,0,0)

    property bool use_orthographic: false
    property bool input_enabled: true
    property bool use_wasd: false

    // Reset targets are scene-owned defaults. The main scene uses a farther
    // orthographic start than the perspective start, while previews can provide
    // their own values if they expose reset controls later.
    property vector3d default_perspective_position: Qt.vector3d(0, 0, 100)
    property vector3d default_orthographic_position: Qt.vector3d(0, 0, 500)
    property quaternion default_camera_rotation: Quaternion.fromEulerAngles(0, 0, 0)

    property real min_camera_distance: 0.01
    property real max_camera_distance: 1000000.0
    property real min_orthographic_magnification: 0.001
    property real max_orthographic_magnification: 1000000.0
    property real default_orthographic_magnification: 100.0

    onEnabledChanged: {
        console.debug("[CameraController] enabled=" + enabled
                      + " input=" + input_enabled
                      + " mouseDown=" + internal.is_mouse_down
                      + " panning=" + internal.is_panning
                      + " mouseAreaPressed=" + cameraMouseArea.pressed)

        if (!enabled && internal.is_mouse_down) {
            mouseReleased(internal.last_pos)
        }
    }

    function clamp_value(value, min_value, max_value) {
        return Math.max(min_value, Math.min(max_value, value))
    }

    function clamp_camera_position(position) {
        var distance = position.length()

        if (distance > max_camera_distance) {
            return position.normalized().times(max_camera_distance)
        }

        return position
    }

    function clamp_orbit_distance(distance) {
        return clamp_value(distance, min_camera_distance, max_camera_distance)
    }

    function clamp_orthographic_magnification(value) {
        return clamp_value(value,
                           min_orthographic_magnification,
                           max_orthographic_magnification)
    }

    // Public entry point used by axis buttons/gizmos. The active navigation
    // mode gets to decide how axis alignment should behave.
    function align_to_axis(axis, invert) {
        internal.current_controller.align_to_axis(axis, invert)
    }

    function align_to_pretty_view() {
        internal.current_controller.align_to_pretty_view()
    }

    function look_at(point) {
        internal.current_controller.look_at(point)
    }

    function reset_view() {
        rotation_target_offset = Qt.vector3d(0, 0, 0)

        perspective_camera.position = clamp_camera_position(
                    default_perspective_position)
        perspective_camera.rotation = default_camera_rotation

        if (orthographic_camera) {
            orthographic_camera.position = clamp_camera_position(
                        default_orthographic_position)
            orthographic_camera.rotation = default_camera_rotation
            orthographic_camera.horizontalMagnification =
                    default_orthographic_magnification
            orthographic_camera.verticalMagnification =
                    default_orthographic_magnification
        }

        internal.current_controller.reset()
    }

    function claim_keyboard_focus() {
        forceActiveFocus()
    }

    function clear_keyboard_input() {
        wasd_control.clear_keyboard_input()
    }

    function as_vector3d(value) {
        return Qt.vector3d(value.x, value.y, value.z)
    }

    function dot(a, b) {
        return a.x * b.x + a.y * b.y + a.z * b.z
    }

    function fit_bounds(bounds_min, bounds_max) {
        var min_point = as_vector3d(bounds_min)
        var max_point = as_vector3d(bounds_max)
        var center = min_point.plus(max_point).times(0.5)
        var extent = max_point.minus(min_point)

        if (extent.x < 0 || extent.y < 0 || extent.z < 0)
            return

        var cam = active_camera
        var forward = as_vector3d(cam.forward).normalized()
        var right = as_vector3d(cam.right).normalized()
        var up = as_vector3d(cam.up).normalized()

        var half_width = 0.0
        var half_height = 0.0
        var half_depth = 0.0

        for (var ix = 0; ix < 2; ++ix) {
            for (var iy = 0; iy < 2; ++iy) {
                for (var iz = 0; iz < 2; ++iz) {
                    var corner = Qt.vector3d(
                                ix ? max_point.x : min_point.x,
                                iy ? max_point.y : min_point.y,
                                iz ? max_point.z : min_point.z)
                    var rel = corner.minus(center)

                    half_width = Math.max(half_width, Math.abs(dot(rel, right)))
                    half_height = Math.max(half_height, Math.abs(dot(rel, up)))
                    half_depth = Math.max(half_depth, Math.abs(dot(rel, forward)))
                }
            }
        }

        var padding = 1.2
        half_width = Math.max(half_width * padding, 1.0)
        half_height = Math.max(half_height * padding, 1.0)
        half_depth = Math.max(half_depth * padding, 1.0)

        var base = rotation_target ? rotation_target.scenePosition : Qt.vector3d(0, 0, 0)
        rotation_target_offset = center.minus(base)

        if (use_orthographic && orthographic_camera) {
            var aspect = Math.max(width, 1) / Math.max(height, 1)
            var horizontal_mag = half_width * 2.0
            var vertical_mag = half_height * 2.0

            if (horizontal_mag / Math.max(vertical_mag, 0.000001) < aspect)
                horizontal_mag = vertical_mag * aspect
            else
                vertical_mag = horizontal_mag / aspect

            orthographic_camera.horizontalMagnification =
                    clamp_orthographic_magnification(horizontal_mag)
            orthographic_camera.verticalMagnification =
                    clamp_orthographic_magnification(vertical_mag)

            var ortho_distance = Math.max(cam.position.minus(center).length(),
                                          half_depth * 4.0,
                                          1.0)
            orthographic_camera.position =
                    clamp_camera_position(center.minus(forward.times(ortho_distance)))
            orthographic_camera.clipFar =
                    Math.max(orthographic_camera.clipFar, ortho_distance + half_depth * 4.0)
        } else {
            var vertical_fov_rad =
                    perspective_camera.fieldOfView * Math.PI / 180.0
            var aspect_ratio = Math.max(width, 1) / Math.max(height, 1)
            var horizontal_fov_rad =
                    2.0 * Math.atan(Math.tan(vertical_fov_rad / 2.0) * aspect_ratio)
            var distance_for_height = half_height / Math.tan(vertical_fov_rad / 2.0)
            var distance_for_width = half_width / Math.tan(horizontal_fov_rad / 2.0)
            var distance = clamp_orbit_distance(
                        Math.max(distance_for_height, distance_for_width) + half_depth)

            perspective_camera.position =
                    clamp_camera_position(center.minus(forward.times(distance)))
            perspective_camera.clipFar =
                    Math.max(perspective_camera.clipFar, distance + half_depth * 4.0)
        }

        internal.current_controller.reset()
    }

    // Converts an enum axis plus an inversion flag into the world-space
    // direction the camera should look from/to for alignment commands.
    function build_align_vector(axis, invert) {
        var axis_setup = Qt.vector3d(0,0,0)

        switch (axis) {
        case CameraController.X:
            axis_setup = Qt.vector3d(1,0,0)
            break;
        case CameraController.Y:
            axis_setup = Qt.vector3d(0,0,1)
            break;
        case CameraController.Z:
            axis_setup = Qt.vector3d(0,1,0)
            break;
        default:
            return
        }

        if (invert) {
            axis_setup = axis_setup.times(-1.0)
        }

        return axis_setup
    }

    function build_pretty_view_vector() {
        return build_align_vector(CameraController.X, false)
            .plus(build_align_vector(CameraController.Y, false))
            .plus(build_align_vector(CameraController.Z, false))
            .normalized()
    }

    onUse_orthographicChanged: {
        if (use_orthographic) {
            use_wasd = false
        }
    }

    onUse_wasdChanged: {
        if (!use_wasd) {
            clear_keyboard_input()
        }
    }

    onInput_enabledChanged: {
        console.debug("[CameraController] input_enabled=" + input_enabled
                      + " enabled=" + enabled
                      + " mouseDown=" + internal.is_mouse_down
                      + " panning=" + internal.is_panning
                      + " mouseAreaPressed=" + cameraMouseArea.pressed)

        if (!input_enabled) {
            clear_keyboard_input()
        }
    }

    // The controller needs focus for keyboard-driven WASD movement. Pointer
    // press below calls forceActiveFocus() so clicking the viewport arms this.
    focus: true

    onActiveFocusChanged: {
        if (!activeFocus) {
            clear_keyboard_input()
        }
    }

    Keys.onPressed: (event)=> { if (input_enabled && !event.isAutoRepeat) handleKeyPress(event) }
    Keys.onReleased: (event)=> { if (input_enabled && !event.isAutoRepeat) handleKeyRelease(event) }

    function handleKeyPress(event) {
        internal.current_controller.handleKeyPress(event)
    }

    function handleKeyRelease(event) {
        internal.current_controller.handleKeyRelease(event)
    }

    MouseArea {
        id: cameraMouseArea
        anchors.fill: parent
        acceptedButtons: Qt.LeftButton
        enabled: root.enabled && root.input_enabled

        onPressed: (mouse) => {
            var pan = Boolean(mouse.modifiers & Qt.ShiftModifier) && !root.use_wasd
            root.mousePressed(Qt.vector2d(mouse.x, mouse.y), pan)
        }

        onPositionChanged: (mouse) => {
            root.mouseMoved(Qt.vector2d(mouse.x, mouse.y))
        }

        onReleased: (mouse) => {
            root.mouseReleased(Qt.vector2d(mouse.x, mouse.y))
        }

        onCanceled: {
            root.mouseReleased(internal.last_pos)
        }
    }

    // ============================================================
    // INTERNAL STATE
    // ============================================================

    QtObject {
        id: internal

        // Needed to guard moves against press and release ordering. Pointer
        // events can arrive around active-state transitions, so we keep
        // explicit state instead of assuming every move belongs to an active
        // drag.
        property bool is_mouse_down: false

        // Latched at mouse press. This prevents mid-drag modifier changes from
        // switching an orbit gesture into a pan gesture halfway through.
        property bool is_panning: false

        // last_pos and mouse_delta_pos are screen-space coordinates in pixels.
        // Individual controllers decide how to convert pixels into angles or
        // world-space movement.
        property vector2d last_pos: Qt.vector2d(0, 0)
        property vector2d mouse_delta_pos: Qt.vector2d(0, 0)

        // All shared input plumbing routes to current_controller, except the
        // Shift-pan branch which is orbit-only.
        property QtObject current_controller : root.use_wasd
                                               ? wasd_control
                                               :
                                                 orbit_control
    }



    function mousePressed(coord, pan) {
        console.debug("[CameraController] mousePressed coord=" + coord
                      + " pan=" + pan
                      + " use_wasd=" + root.use_wasd)
        forceActiveFocus()

        // The gesture mode is captured here and preserved until release. That
        // keeps the drag stable if the user presses/releases Shift while moving.
        internal.is_mouse_down = true
        internal.is_panning = pan
        internal.mouse_delta_pos = Qt.vector2d(0,0)
        internal.last_pos = coord

        if (!root.use_wasd) {
            // Orbit math is stateful. Rebuild yaw/pitch from the camera at the
            // beginning of each gesture so external camera changes, axis-align
            // animations, or prior pans are reflected before the next drag.
            orbit_control.initialize_from_camera()
        }
    }

    function mouseReleased(coord) {
        console.debug("[CameraController] mouseReleased coord=" + coord
                      + " wasMouseDown=" + internal.is_mouse_down
                      + " wasPanning=" + internal.is_panning)
        internal.is_mouse_down = false
        internal.is_panning = false
        internal.mouse_delta_pos = Qt.vector2d(0,0)
        internal.last_pos = coord
    }

    function mouseMoved(coord) {
        if (!internal.is_mouse_down) {
            console.debug("[CameraController] mouseMoved ignored coord=" + coord
                          + " mouseDown=false")
            return
        }

        internal.mouse_delta_pos = coord.minus(internal.last_pos)
        internal.last_pos = coord

        console.debug("[CameraController] mouseMoved coord=" + coord
                      + " delta=" + internal.mouse_delta_pos
                      + " panning=" + internal.is_panning
                      + " use_wasd=" + root.use_wasd)

        // Pan is handled outside current_controller because it is not a general
        // navigation mode. It is an orbit gesture modifier that translates the
        // orbit center and camera together.
        if (internal.is_panning && !root.use_wasd)
            orbit_control.apply_pan_delta(coord)
        else
            internal.current_controller.apply_mouse_delta(coord)

        internal.mouse_delta_pos = Qt.vector2d(0,0)
    }

    WheelHandler {
        id: wheel_handler

        orientation: Qt.Vertical
        target: null
        enabled: root.enabled && root.input_enabled

        acceptedDevices: PointerDevice.Mouse | PointerDevice.TouchPad

        onWheel: function(event) {
            if (root.use_wasd) {
                // In fly mode, scroll changes movement speed instead of moving
                // the camera immediately. This mirrors common 3D editor camera
                // controls and keeps fine/coarse navigation accessible.
                wasd_control.speed_multiplier *=
                        Math.exp(wasd_control.scroll_factor * event.angleDelta.y)

                wasd_control.speed_multiplier = Math.max(
                            Math.min(wasd_control.speed_multiplier, 100000),
                            0.000001)
            } else {
                orbit_control.apply_wheel_delta(event.angleDelta.y)
            }
        }
    }

    FrameAnimation {
        running: true
        onTriggered: {
            // Continuous controllers advance here. Orbit mostly reacts to
            // discrete pointer/wheel events, but axis alignment uses this frame
            // callback for smooth animation.
            internal.current_controller.processInput(frameTime)
        }
    }

    // Shared animation object for WASD axis alignment. Orbit alignment uses its
    // own scalar interpolation because the orbit state is yaw/pitch/distance
    // rather than raw camera position.
    ParallelAnimation {
        id: wasd_camera_animation

        Vector3dAnimation {
            id: camera_move_animation
            duration: 250

            target: root.active_camera
            property: "position"

            easing: Easing.InOutCubic

        }

        QuaternionAnimation {
            id: camera_rotation_animation
            duration: 250

            target: root.active_camera
            property: "rotation"

            easing: Easing.InOutCubic
        }

        onFinished: {
            wasd_control.is_animating = false
        }
    }

    CameraWasdController {
        id: wasd_control
        host: root
        inputState: internal
        cameraMoveAnimation: camera_move_animation
        cameraRotationAnimation: camera_rotation_animation
        cameraAnimation: wasd_camera_animation
    }

    CameraOrbitController {
        id: orbit_control
        host: root
        inputState: internal
    }
}
