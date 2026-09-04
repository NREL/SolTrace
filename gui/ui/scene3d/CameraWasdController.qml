import QtQuick
import QtQuick3D

import SolTrace

QtObject {
    id: root

    required property var host
    required property var inputState
    required property var cameraMoveAnimation
    required property var cameraRotationAnimation
    required property var cameraAnimation

    // WASD mode is camera-relative free flight:
    //   W/S move forward/back along the camera forward vector.
    //   A/D strafe along the camera right vector.
    //   Q/E move down/up in world Y.
    //   Shift increases speed.
    //
    // This mode intentionally ignores rotation_target.
    property real sensitivity: 0.2
    property real walk_speed: App.view.sim.fps_walk_speed
    property real run_speed: walk_speed + 10

    // Multiplicative scroll scale. The constant produces gentle exponential
    // changes from wheel deltas without needing frame-rate dependent input.
    property real scroll_factor: 0.04879016 / 4.0

    // Applied only when no movement keys are active, so the camera eases to
    // a stop instead of snapping instantly.
    property real friction: 1.5
    property real mouse_sensitivity_deg: 0.3

    property bool is_animating: false
    property bool initialized: false
    property real speed_multiplier: 1.0
    property vector3d velocity: Qt.vector3d(0, 0, 0)

    property vector3d kb_state: Qt.vector3d(0, 0, 0)
    property bool kb_shift: false

    function negate(vector) {
        return Qt.vector3d(-vector.x, -vector.y, -vector.z)
    }

    function processInput(frameDelta) {
        if (is_animating) {
            return
        }

        var is_moving = !kb_state.fuzzyEquals(Qt.vector3d(0, 0, 0))

        if (is_moving) {
            let speed = (kb_shift ? run_speed : walk_speed) * speed_multiplier

            velocity = kb_state.normalized().times(speed)
        } else {
            velocity = velocity.times(1.0 / friction)

            if (velocity.length() < 1e-6) {
                velocity = Qt.vector3d(0, 0, 0)
            }
        }

        is_moving = !velocity.fuzzyEquals(Qt.vector3d(0, 0, 0))

        if (is_moving) {
            // Qt exposes camera basis vectors as QVector3D-like values. We
            // copy them into QML vector3d values before using arithmetic
            // helpers such as times()/plus().
            let forward = host.active_camera.forward
            let right = host.active_camera.right

            forward = Qt.vector3d(forward.x, forward.y, forward.z)
            right = Qt.vector3d(right.x, right.y, right.z)

            let mlt = velocity.times(frameDelta)

            var xp = right.times(mlt.x)
            var yp = Qt.vector3d(0, 1, 0).times(mlt.y)
            var zp = forward.times(mlt.z)

            let delta = xp.plus(yp).plus(zp)

            let current_pos = host.active_camera.position

            let new_pos = Qt.vector3d(current_pos.x,
                                      current_pos.y,
                                      current_pos.z).plus(delta)

            host.active_camera.position = host.clamp_camera_position(new_pos)
        }
    }

    function apply_mouse_delta(coord) {
        // Free-look changes camera Euler angles directly. This is separate
        // from orbit mode, which recomputes camera position from yaw/pitch
        // around a target point.
        var rotationVector = host.active_camera.eulerRotation
        let pitch = rotationVector.x
        let yaw = rotationVector.y

        pitch = (pitch - inputState.mouse_delta_pos.y * mouse_sensitivity_deg * sensitivity)
        pitch = Math.max(Math.min(pitch, 89), -89)

        yaw -= inputState.mouse_delta_pos.x * mouse_sensitivity_deg * sensitivity

        host.active_camera.eulerRotation.x = pitch
        host.active_camera.eulerRotation.y = yaw
    }

    function rotation_from_forward(direction) {
        var forward = direction.normalized()
        var horizontal = Math.sqrt(forward.x * forward.x
                                   + forward.z * forward.z)
        var yaw = Math.atan2(-forward.x, -forward.z) * 180.0 / Math.PI
        var pitch = Math.atan2(forward.y, horizontal) * 180.0 / Math.PI

        pitch = Math.max(Math.min(pitch, 89), -89)

        return Quaternion.fromEulerAngles(pitch, yaw, 0)
    }

    function align_to_axis(axis, invert) {
        var axis_setup = host.build_align_vector(axis, invert)

        // WASD alignment rotates in place. Position is animated from/to the
        // same value so the shared ParallelAnimation can drive both
        // position and rotation targets without a special case.
        var rotation = rotation_from_forward(axis_setup)

        var camera_position = host.clamp_camera_position(
                    host.active_camera.position)

        cameraMoveAnimation.from = camera_position
        cameraMoveAnimation.to = camera_position

        cameraRotationAnimation.from = host.active_camera.rotation
        cameraRotationAnimation.to = rotation

        is_animating = true
        cameraAnimation.start()
    }

    function align_to_pretty_view() {
        var rotation = rotation_from_forward(
                    host.build_pretty_view_vector().times(-1.0))

        var camera_position = host.clamp_camera_position(
                    host.active_camera.position)

        cameraMoveAnimation.from = camera_position
        cameraMoveAnimation.to = camera_position

        cameraRotationAnimation.from = host.active_camera.rotation
        cameraRotationAnimation.to = rotation

        is_animating = true
        cameraAnimation.start()
    }

    function look_at(point) {
        var cam = host.active_camera
        var camera_position = host.clamp_camera_position(cam.position)
        var target = Qt.vector3d(point.x, point.y, point.z)
        var offset = target.minus(camera_position)

        if (offset.length() < 0.000001)
            return

        cameraMoveAnimation.from = camera_position
        cameraMoveAnimation.to = camera_position

        cameraRotationAnimation.from = cam.rotation
        cameraRotationAnimation.to = rotation_from_forward(offset)

        is_animating = true
        cameraAnimation.start()
    }

    function clear_keyboard_input() {
        velocity = Qt.vector3d(0, 0, 0)
        kb_state = Qt.vector3d(0, 0, 0)
        kb_shift = false
    }

    function reset() {
        is_animating = false
        clear_keyboard_input()
        speed_multiplier = 1.0
        cameraAnimation.stop()
    }

    function handleKeyPress(event) {
        switch (event.key) {
        case Qt.Key_W:
            kb_state.z = 1.0
            break
        case Qt.Key_S:
            kb_state.z = -1.0
            break
        case Qt.Key_A:
            kb_state.x = -1.0
            break
        case Qt.Key_D:
            kb_state.x = 1.0
            break
        case Qt.Key_Q:
            kb_state.y = -1.0
            break
        case Qt.Key_E:
            kb_state.y = 1.0
            break
        case Qt.Key_Shift:
            kb_shift = true
            break
        }
    }

    function handleKeyRelease(event) {
        switch (event.key) {
        case Qt.Key_W:
        case Qt.Key_S:
            kb_state.z = 0.0
            break
        case Qt.Key_A:
        case Qt.Key_D:
            kb_state.x = 0.0
            break
        case Qt.Key_Q:
        case Qt.Key_E:
            kb_state.y = 0.0
            break
        case Qt.Key_Shift:
            kb_shift = false
            break
        }
    }
}
