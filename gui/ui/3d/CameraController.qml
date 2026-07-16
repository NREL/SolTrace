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

    // Public entry point used by axis buttons/gizmos. The active navigation
    // mode gets to decide how axis alignment should behave.
    function align_to_axis(axis, invert) {
        internal.current_controller.align_to_axis(axis, invert)
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
            axis_setup = Qt.vector3d(0,1,0)
            break;
        case CameraController.Z:
            axis_setup = Qt.vector3d(0,0,1)
            break;
        default:
            return
        }

        if (invert) {
            axis_setup = axis_setup.times(-1.0)
        }

        return axis_setup
    }

    onUse_orthographicChanged: {
        if (use_orthographic) {
            use_wasd = false
        }
    }

    // The controller needs focus for keyboard-driven WASD movement. Pointer
    // press below calls forceActiveFocus() so clicking the viewport arms this.
    focus: true

    Keys.onPressed: (event)=> { if (input_enabled && !event.isAutoRepeat) handleKeyPress(event) }
    Keys.onReleased: (event)=> { if (input_enabled && !event.isAutoRepeat) handleKeyRelease(event) }

    function handleKeyPress(event) {
        internal.current_controller.handleKeyPress(event)
    }

    function handleKeyRelease(event) {
        internal.current_controller.handleKeyRelease(event)
    }

    DragHandler {
        id: dragHandler
        target: null
        enabled: root.input_enabled

        // Plain drag means "rotate/look". Keeping this handler NoModifier lets
        // the Shift-specific handler below own panning without mode checks in
        // the gesture recognizer itself.
        acceptedModifiers: Qt.NoModifier

        onCentroidChanged: {
            root.mouseMoved(Qt.vector2d(centroid.position.x, centroid.position.y), false);
        }

        onActiveChanged: {
            if (active)
                root.mousePressed(Qt.vector2d(centroid.position.x, centroid.position.y), false);
            else
                root.mouseReleased(Qt.vector2d(centroid.position.x, centroid.position.y));
        }
    }

    DragHandler {
        id: panDragHandler
        target: null
        enabled: root.input_enabled && !root.use_wasd

        // Orbit panning is intentionally tied to Shift-drag. In WASD mode Shift
        // already means "run", so pan is disabled there to avoid overloading
        // the same gesture.
        acceptedModifiers: Qt.ShiftModifier

        onCentroidChanged: {
            root.mouseMoved(Qt.vector2d(centroid.position.x, centroid.position.y), true);
        }

        onActiveChanged: {
            if (active)
                root.mousePressed(Qt.vector2d(centroid.position.x, centroid.position.y), true);
            else
                root.mouseReleased(Qt.vector2d(centroid.position.x, centroid.position.y));
        }
    }

    // ============================================================
    // INTERNAL STATE
    // ============================================================

    QtObject {
        id: internal

        // Needed to guard moves against press and release ordering. DragHandler
        // can emit centroid changes around active-state transitions, so we keep
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
        internal.is_mouse_down = false
        internal.is_panning = false
        internal.mouse_delta_pos = Qt.vector2d(0,0)
        internal.last_pos = coord
    }

    function mouseMoved(coord) {
        if (!internal.is_mouse_down) return

        internal.mouse_delta_pos = coord.minus(internal.last_pos)
        internal.last_pos = coord

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
        enabled: root.input_enabled

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

    // ============================================================
    // WASD CONTROLLER
    // ============================================================

    QtObject {
        id: wasd_control

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
        property vector3d velocity: Qt.vector3d(0,0,0)

        property vector3d kb_state: Qt.vector3d(0,0,0)
        property bool kb_shift: false


        function negate(vector) {
            return Qt.vector3d(-vector.x, -vector.y, -vector.z)
        }

        function processInput(frameDelta) {
            if (is_animating) {
                return;
            }

            var is_moving = !kb_state.fuzzyEquals(Qt.vector3d(0,0,0))

            if (is_moving) {
                let speed = (kb_shift ? run_speed : walk_speed) * speed_multiplier

                velocity = kb_state.normalized().times(speed)
            } else {
                velocity = velocity.times(1.0 / friction)

                if (velocity.length() < 1e-6) {
                    velocity = Qt.vector3d(0,0,0)
                }
            }

            is_moving = !velocity.fuzzyEquals(Qt.vector3d(0,0,0))

            if (is_moving) {
                // Qt exposes camera basis vectors as QVector3D-like values. We
                // copy them into QML vector3d values before using arithmetic
                // helpers such as times()/plus().
                let forward = root.active_camera.forward
                let right = root.active_camera.right

                forward = Qt.vector3d(forward.x, forward.y, forward.z)
                right = Qt.vector3d(right.x, right.y, right.z)

                let mlt = velocity.times(frameDelta)

                var xp = right.times(mlt.x)
                var yp = Qt.vector3d(0,1,0).times(mlt.y)
                var zp = forward.times(mlt.z)

                let delta = xp.plus(yp).plus(zp)

                let current_pos = root.active_camera.position

                let new_pos = Qt.vector3d(current_pos.x,
                                          current_pos.y,
                                          current_pos.z).plus(delta)

                root.active_camera.position = root.clamp_camera_position(new_pos)
            }
        }

        function apply_mouse_delta(coord) {
            // Free-look changes camera Euler angles directly. This is separate
            // from orbit mode, which recomputes camera position from yaw/pitch
            // around a target point.
            var rotationVector = root.active_camera.eulerRotation;
            let pitch = rotationVector.x
            let yaw = rotationVector.y

            pitch = (pitch - internal.mouse_delta_pos.y * mouse_sensitivity_deg * sensitivity)
            pitch = Math.max(Math.min(pitch, 89), - 89)

            yaw -= internal.mouse_delta_pos.x * mouse_sensitivity_deg * sensitivity

            root.active_camera.eulerRotation.x = pitch
            root.active_camera.eulerRotation.y = yaw
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
            var axis_setup = root.build_align_vector(axis, invert)

            // WASD alignment rotates in place. Position is animated from/to the
            // same value so the shared ParallelAnimation can drive both
            // position and rotation targets without a special case.
            var rotation = rotation_from_forward(axis_setup)

            var camera_position = root.clamp_camera_position(
                        root.active_camera.position)

            camera_move_animation.from = camera_position
            camera_move_animation.to = camera_position

            camera_rotation_animation.from = root.active_camera.rotation
            camera_rotation_animation.to = rotation

            is_animating = true
            wasd_camera_animation.start()
        }

        function look_at(point) {
            var cam = root.active_camera
            var camera_position = root.clamp_camera_position(cam.position)
            var target = Qt.vector3d(point.x, point.y, point.z)
            var offset = target.minus(camera_position)

            if (offset.length() < 0.000001)
                return

            camera_move_animation.from = camera_position
            camera_move_animation.to = camera_position

            camera_rotation_animation.from = cam.rotation
            camera_rotation_animation.to = rotation_from_forward(offset)

            is_animating = true
            wasd_camera_animation.start()
        }

        function reset() {
            is_animating = false
            velocity = Qt.vector3d(0, 0, 0)
            kb_state = Qt.vector3d(0, 0, 0)
            kb_shift = false
            speed_multiplier = 1.0
            wasd_camera_animation.stop()
        }

        function handleKeyPress(event) {
            switch (event.key) {
            case Qt.Key_W:
                kb_state.z = 1.0
                break;
            case Qt.Key_S:
                kb_state.z = -1.0
                break;
            case Qt.Key_A:
                kb_state.x = -1.0
                break;
            case Qt.Key_D:
                kb_state.x = 1.0
                break;
            case Qt.Key_Q:
                kb_state.y = -1.0
                break;
            case Qt.Key_E:
                kb_state.y = 1.0
                break;
            case Qt.Key_Shift:
                kb_shift = true
                break;
            }
        }

        function handleKeyRelease(event) {
            switch (event.key) {
            case Qt.Key_W:
            case Qt.Key_S:
                kb_state.z = 0.0
                break;
            case Qt.Key_A:
            case Qt.Key_D:
                kb_state.x = 0.0
                break;
            case Qt.Key_Q:
            case Qt.Key_E:
                kb_state.y = 0.0
                break;
            case Qt.Key_Shift:
                kb_shift = false
                break;
            }
        }
    }

    // ============================================================
    // ORBIT CONTROLLER
    // ============================================================

    QtObject {
        id: orbit_control

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
            let base = root.rotation_target
                ? root.rotation_target.scenePosition
                : Qt.vector3d(0, 0, 0)

            return base.plus(root.rotation_target_offset)
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

            apply_orbit_transform(root.clamp_orbit_distance(distance))

            if (t >= 1.0) {
                is_animating = false
                yaw_deg = animation_to_yaw_deg
                pitch_deg = animation_to_pitch_deg
                apply_orbit_transform(root.clamp_orbit_distance(
                                          animation_to_distance))
            }
        }

        function initialize_from_camera() {
            var cam = root.active_camera
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
                        Qt.vector3d(0,1,0),
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
            var cam = root.active_camera
            var target = rotation_point

            var q = orbit_rotation()

            // In this convention the camera sits on the local +Z axis at the
            // requested distance, then the orbit rotation moves that offset into
            // world space. The camera rotation is the same quaternion so it
            // looks back toward the target.
            var orbit_distance = root.clamp_orbit_distance(distance)
            var local_orbit_offset = Qt.vector3d(0, 0, orbit_distance)
            var new_offset = q.times(local_orbit_offset)

            cam.position = root.clamp_camera_position(target.plus(new_offset))
            cam.rotation = q
        }


        function apply_mouse_delta(coord) {
            var cam = root.active_camera
            var target = rotation_point

            // Preserve the current radius while changing angular state.
            var offset = cam.position.minus(target)
            var distance = offset.length()

            if (distance < 0.000001)
                distance = 0.000001

            yaw_deg -= internal.mouse_delta_pos.x *
                    mouse_sensitivity_deg *
                    sensitivity

            pitch_deg += internal.mouse_delta_pos.y *
                    mouse_sensitivity_deg *
                    sensitivity

            pitch_deg = clamp(pitch_deg, min_pitch_deg, max_pitch_deg)

            apply_orbit_transform(distance)
        }

        function apply_pan_delta(coord) {
            var cam = root.active_camera
            var target = rotation_point

            // Pan is a screen-space translation mapped to camera-right and
            // camera-up. The orbit center and the camera move together, which
            // keeps the apparent view direction unchanged while sliding the
            // scene under the cursor.
            var offset = cam.position.minus(target)
            var distance = offset.length()

            if (distance < 0.000001)
                distance = 1.0

            var view_size = Math.max(1, Math.min(root.width, root.height))

            // Scale by distance so panning feels similar regardless of zoom.
            // Near the target, small pixel moves become fine world-space moves;
            // far away, the same gesture covers more world space.
            var world_units_per_pixel = distance / view_size * pan_sensitivity

            // Dragging right should move the scene right on screen, which means
            // the camera/target move left in world camera-right space. Dragging
            // down should move the scene down, so camera/target move up.
            var right = Qt.vector3d(cam.right.x, cam.right.y, cam.right.z).normalized()
            var up = Qt.vector3d(cam.up.x, cam.up.y, cam.up.z).normalized()

            var delta = right.times(-internal.mouse_delta_pos.x * world_units_per_pixel)
                .plus(up.times(internal.mouse_delta_pos.y * world_units_per_pixel))

            var requested_position = cam.position.plus(delta)
            var clamped_position = root.clamp_camera_position(requested_position)
            var applied_delta = clamped_position.minus(cam.position)

            root.rotation_target_offset =
                    root.rotation_target_offset.plus(applied_delta)
            cam.position = clamped_position
        }


        function apply_wheel_delta(delta_y) {
            var cam = root.active_camera
            var target = rotation_point

            // Refresh yaw/pitch before zooming so wheel input after external
            // camera movement keeps orbit state synchronized.
            initialize_from_camera()

            var offset = cam.position.minus(target)
            var distance = offset.length()

            if (distance < 0.000001)
                return

            var zoom_scale = Math.exp(-delta_y * zoom_factor)
            var new_distance = root.clamp_orbit_distance(distance * zoom_scale)

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

        function align_to_axis(axis, invert) {
            var cam = root.active_camera
            var target = rotation_point

            // Keep the current distance and animate only the angular state.
            initialize_from_camera()

            var current_offset = cam.position.minus(target)
            var current_distance = current_offset.length()

            if (current_distance < 0.000001)
                current_distance = 1.0

            current_distance = root.clamp_orbit_distance(current_distance)

            var axis_offset = build_align_vector(axis, invert)
            var desired_offset = axis_offset.times(current_distance)

            var angles = yaw_pitch_from_offset(desired_offset)

            var target_pitch = clamp(angles.pitch,
                                     min_pitch_deg,
                                     max_pitch_deg)

            start_orbit_animation(angles.yaw,
                                  target_pitch,
                                  current_distance)
        }

        function look_at(point) {
            var cam = root.active_camera
            var target = Qt.vector3d(point.x, point.y, point.z)
            var base = root.rotation_target
                    ? root.rotation_target.scenePosition
                    : Qt.vector3d(0, 0, 0)

            root.rotation_target_offset = target.minus(base)

            var offset = cam.position.minus(rotation_point)
            var distance = offset.length()

            if (distance < 0.000001)
                return

            distance = root.clamp_orbit_distance(distance)

            is_animating = false
            initialize_from_camera()
            apply_orbit_transform(distance)
        }

        function start_orbit_animation(to_yaw_deg, to_pitch_deg, to_distance) {
            var cam = root.active_camera
            var target = rotation_point

            // Capture the current orbit state as the animation start point.
            initialize_from_camera()

            var offset = cam.position.minus(target)
            var distance = offset.length()

            if (distance < 0.000001)
                distance = 1.0

            distance = root.clamp_orbit_distance(distance)

            animation_from_yaw_deg = yaw_deg
            animation_from_pitch_deg = pitch_deg
            animation_from_distance = distance

            animation_to_yaw_deg = to_yaw_deg
            animation_to_pitch_deg = to_pitch_deg
            animation_to_distance = root.clamp_orbit_distance(to_distance)

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
            // by panDragHandler as a pointer modifier rather than here.
        }

        function handleKeyRelease(event) {
            // Orbit mode has no key-release state to clear.
        }
    }
}
