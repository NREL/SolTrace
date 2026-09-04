# QML Source Organization

The QML tree is organized by dependency direction:

- `app/`: application entry points and global app-level objects.
- `shell/`: persistent application chrome, workflow navigation, top bar, bottom bar, and panel hosts.
- `features/`: workflow/domain UI. Feature code may depend on controls, fields, widgets, shell state, `App`, and `AppData`.
- `scene3d/`: Qt Quick 3D scene, camera, picking, and gizmo UI.
- `controls/`: design-system primitives such as buttons, combo boxes, form rows, popups, and panel containers. These should stay broadly reusable.
- `fields/`: reusable field composites such as location, date/time, slider, and color pickers.
- `widgets/`: reusable app widgets that may know about application models or workflows.

Prefer placing new QML near the feature that owns it. Move code into `controls/`,
`fields/`, or `widgets/` only when it is reused or clearly generic.
