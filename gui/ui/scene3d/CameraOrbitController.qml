import QtQuick
import QtQuick3D

QtObject {
    id: root

    required property var host
    required property var inputState

    // Orbit mode stores camera orientation as yaw/pitch plus a distance to
    // rotation_point. Every orbit/zoom/axis operation reconstructs camera
    // position from that state. This avoids accumulating small transform
    // errors from repeated relative rotations.
    property real sensitivity: 0.2
    property real mouse_sensitivity_deg: 0.6

    // yaw_deg rotates around world Y. pitch_deg tilts above/below the
    // target; positive pitch means the camera is above the target.
    property real yaw_deg: 0
    property real pitch_deg: 0

    // Avoid the singularity where the camera is exactly above/below the
    // target and the right axis becomes unstable.
    property real min_pitch_deg: -89
    property real max_pitch_deg: 89

    // Wheel zoom is exponential, so it feels similar at small and large
    // distances. Pan sensitivity is applied after converting pixels to
    // world units.
    property real zoom_factor: 0.0005
    property real pan_sensitivity: 1.0

    property bool is_animating: false

    property real animation_duration: 0.25
    property real animation_time: 0.0

    property real animation_from_yaw_deg: 0
    property real animation_from_pitch_deg: 0
    property real animation_from_distance: 1

    property real animation_to_yaw_deg: 0
    property real animation_to_pitch_deg: 0
    property real animation_to_distance: 1

    // The effective center of orbit. If a scene target exists, start from
    // its scene position; otherwise, use world origin. Panning contributes
    // the accumulated offset without mutating the target node.
    property vector3d rotation_point: {
        let base = host.rotation_target
            ? host.rotation_target.scenePosition
            : Qt.vector3d(0, 0, 0)

        return base.plus(host.rotation_target_offset)
    }

    function lerp(a, b, t) {
        return a + (b - a) * t
    }

    function shortest_angle_delta_deg(from_deg, to_deg) {
        var delta = to_deg - from_deg

        while (delta > 180.0)
            delta -= 360.0

        while (delta < -180.0)
            delta += 360.0

        return delta
    }

    function lerp_angle_deg(from_deg, to_deg, t) {
        return from_deg + shortest_angle_delta_deg(from_deg, to_deg) * t
    }

    function deg_to_rad(deg) {
        return deg * Math.PI / 180.0
    }

    function clamp(value, min_value, max_value) {
        return Math.max(min_value, Math.min(max_value, value))
    }

    function processInput(frameDelta) {
        if (!is_animating)
            return

        // Axis alignment is animated in orbit-space coordinates so the
        // camera stays on the same orbit radius while yaw/pitch ease to the
        // requested axis direction.
        animation_time += frameDelta

        var t = animation_time / animation_duration
        t = clamp(t, 0.0, 1.0)

        // Smoothstep easing.
        var eased = t * t * (3.0 - 2.0 * t)

        yaw_deg = lerp_angle_deg(animation_from_yaw_deg,
                                 animation_to_yaw_deg,
                                 eased)

        pitch_deg = lerp(animation_from_pitch_deg,
                         animation_to_pitch_deg,
                         eased)

        var distance = lerp(animation_from_distance,
                            animation_to_distance,
                            eased)

        apply_orbit_transform(host.clamp_orbit_distance(distance))

        if (t >= 1.0) {
            is_animating = false
            yaw_deg = animation_to_yaw_deg
            pitch_deg = animation_to_pitch_deg
            apply_orbit_transform(host.clamp_orbit_distance(
                                      animation_to_distance))
        }
    }

    function initialize_from_camera() {
        var cam = host.active_camera
        var target = rotation_point

        // Convert current camera position into orbit state. This is called
        // before user gestures and animations because other code may have
        // moved the camera since the last orbit update.
        var offset = cam.position.minus(target)

        var horizontal_distance = Math.sqrt(
                    offset.x * offset.x +
                    offset.z * offset.z)

        yaw_deg = Math.atan2(offset.x, offset.z) * 180.0 / Math.PI
        pitch_deg = Math.atan2(offset.y, horizontal_distance) * 180.0 / Math.PI

        pitch_deg = clamp(pitch_deg, min_pitch_deg, max_pitch_deg)
    }

    function orbit_rotation() {
        // Build the camera rotation that corresponds to the stored orbit
        // yaw/pitch. Yaw rotates around world up first. The pitch axis then
        // becomes the yawed camera-right vector.
        var yaw_q = Quaternion.fromAxisAndAngle(
                    Qt.vector3d(0, 1, 0),
                    yaw_deg)

        var right_axis = yaw_q.times(Qt.vector3d(1, 0, 0)).normalized()

        // Important: negate pitch here.
        // initialize_from_camera() reads positive pitch as camera above target.
        // Qt rotation around +X needs negative angle to recreate that offset.
        var pitch_q = Quaternion.fromAxisAndAngle(
                    right_axis,
                    -pitch_deg)

        return pitch_q.times(yaw_q).normalized()
    }

    function apply_orbit_transform(distance) {
        var cam = host.active_camera
        var target = rotation_point

        var q = orbit_rotation()

        // In this convention the camera sits on the local +Z axis at the
        // requested distance, then the orbit rotation moves that offset into
        // world space. The camera rotation is the same quaternion so it
        // looks back toward the target.
        var orbit_distance = host.clamp_orbit_distance(distance)
        var local_orbit_offset = Qt.vector3d(0, 0, orbit_distance)
        var new_offset = q.times(local_orbit_offset)

        cam.position = host.clamp_camera_position(target.plus(new_offset))
        cam.rotation = q
    }

    function apply_mouse_delta(coord) {
        var cam = host.active_camera
        var target = rotation_point

        // Preserve the current radius while changing angular state.
        var offset = cam.position.minus(target)
        var distance = offset.length()

        if (distance < 0.000001)
            distance = 0.000001

        yaw_deg -= inputState.mouse_delta_pos.x *
                mouse_sensitivity_deg *
                sensitivity

        pitch_deg += inputState.mouse_delta_pos.y *
                mouse_sensitivity_deg *
                sensitivity

        pitch_deg = clamp(pitch_deg, min_pitch_deg, max_pitch_deg)

        apply_orbit_transform(distance)
    }

    function apply_pan_delta(coord) {
        var cam = host.active_camera
        var target = rotation_point

        // Pan is a screen-space translation mapped to camera-right and
        // camera-up. The orbit center and the camera move together, which
        // keeps the apparent view direction unchanged while sliding the
        // scene under the cursor.
        var offset = cam.position.minus(target)
        var distance = offset.length()

        if (distance < 0.000001)
            distance = 1.0

        var view_size = Math.max(1, Math.min(host.width, host.height))

        var visible_size = host.use_orthographic && host.orthographic_camera
                ? Math.max(host.orthographic_camera.horizontalMagnification,
                           host.orthographic_camera.verticalMagnification)
                : distance

        // Scale by visible size so panning feels similar regardless of
        // perspective distance or orthographic zoom.
        var world_units_per_pixel = visible_size / view_size * pan_sensitivity

        // Dragging right should move the scene right on screen, which means
        // the camera/target move left in world camera-right space. Dragging
        // down should move the scene down, so camera/target move up.
        var right = Qt.vector3d(cam.right.x, cam.right.y, cam.right.z).normalized()
        var up = Qt.vector3d(cam.up.x, cam.up.y, cam.up.z).normalized()

        var delta = right.times(-inputState.mouse_delta_pos.x * world_units_per_pixel)
            .plus(up.times(inputState.mouse_delta_pos.y * world_units_per_pixel))

        var requested_position = cam.position.plus(delta)
        var clamped_position = host.clamp_camera_position(requested_position)
        var applied_delta = clamped_position.minus(cam.position)

        host.rotation_target_offset =
                host.rotation_target_offset.plus(applied_delta)
        cam.position = clamped_position
    }

    function apply_wheel_delta(delta_y) {
        var cam = host.active_camera
        var target = rotation_point

        if (host.use_orthographic && host.orthographic_camera) {
            var zoom_scale = Math.exp(-delta_y * zoom_factor)
            var next_magnification = host.clamp_orthographic_magnification(
                        host.orthographic_camera.horizontalMagnification
                        * zoom_scale)

            host.orthographic_camera.horizontalMagnification =
                    next_magnification
            host.orthographic_camera.verticalMagnification =
                    next_magnification
            return
        }

        // Refresh yaw/pitch before zooming so wheel input after external
        // camera movement keeps orbit state synchronized.
        initialize_from_camera()

        var offset = cam.position.minus(target)
        var distance = offset.length()

        if (distance < 0.000001)
            return

        var zoom_scale = Math.exp(-delta_y * zoom_factor)
        var new_distance = host.clamp_orbit_distance(distance * zoom_scale)

        apply_orbit_transform(new_distance)
    }

    function yaw_pitch_from_offset(offset) {
        // Helper for axis alignment: given a desired camera offset from the
        // target, calculate the orbit angles needed to reach it.
        var horizontal_distance = Math.sqrt(
                    offset.x * offset.x +
                    offset.z * offset.z)

        var result = {
            yaw: Math.atan2(offset.x, offset.z) * 180.0 / Math.PI,
            pitch: Math.atan2(offset.y, horizontal_distance) * 180.0 / Math.PI
        }

        return result
    }

    function align_to_offset(offset) {
        var cam = host.active_camera
        var target = rotation_point

        // Keep the current distance and animate only the angular state.
        initialize_from_camera()

        var current_offset = cam.position.minus(target)
        var current_distance = current_offset.length()

        if (current_distance < 0.000001)
            current_distance = 1.0

        current_distance = host.clamp_orbit_distance(current_distance)

        var desired_offset = offset.normalized().times(current_distance)

        var angles = yaw_pitch_from_offset(desired_offset)

        var target_pitch = clamp(angles.pitch,
                                 min_pitch_deg,
                                 max_pitch_deg)

        start_orbit_animation(angles.yaw,
                              target_pitch,
                              current_distance)
    }

    function align_to_axis(axis, invert) {
        align_to_offset(host.build_align_vector(axis, invert))
    }

    function align_to_pretty_view() {
        align_to_offset(host.build_pretty_view_vector())
    }

    function look_at(point) {
        var cam = host.active_camera
        var target = Qt.vector3d(point.x, point.y, point.z)
        var base = host.rotation_target
                ? host.rotation_target.scenePosition
                : Qt.vector3d(0, 0, 0)

        host.rotation_target_offset = target.minus(base)

        var offset = cam.position.minus(rotation_point)
        var distance = offset.length()

        if (distance < 0.000001)
            return

        distance = host.clamp_orbit_distance(distance)

        is_animating = false
        initialize_from_camera()
        apply_orbit_transform(distance)
    }

    function start_orbit_animation(to_yaw_deg, to_pitch_deg, to_distance) {
        var cam = host.active_camera
        var target = rotation_point

        // Capture the current orbit state as the animation start point.
        initialize_from_camera()

        var offset = cam.position.minus(target)
        var distance = offset.length()

        if (distance < 0.000001)
            distance = 1.0

        distance = host.clamp_orbit_distance(distance)

        animation_from_yaw_deg = yaw_deg
        animation_from_pitch_deg = pitch_deg
        animation_from_distance = distance

        animation_to_yaw_deg = to_yaw_deg
        animation_to_pitch_deg = to_pitch_deg
        animation_to_distance = host.clamp_orbit_distance(to_distance)

        animation_time = 0.0
        is_animating = true
    }

    function reset() {
        is_animating = false
        animation_time = 0.0
        initialize_from_camera()
    }

    function handleKeyPress(event) {
        // Orbit mode does not consume keyboard movement. Shift is handled
        // by the camera mouse area as a pointer modifier rather than here.
    }

    function handleKeyRelease(event) {
        // Orbit mode has no key-release state to clear.
    }
}
