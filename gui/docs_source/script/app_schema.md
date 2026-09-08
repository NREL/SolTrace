# App Script API

These globals control application state outside the current database. The database editing API remains available as `db`.

## sim.start(config)

Start a simulation using the current scene and current runner settings, patched by the optional `config` object.

Unspecified keys use the current UI defaults. The function returns `true` when a run was started and `false` when no current scene is available, a simulation is already running, or the selected configuration could not start.

Supported configuration keys:

- `runner`: runner name or alias. Accepted aliases include `cpu`, `legacy`, `native`, `embree`, `gpu`, and `optix`.
- `runner_index`: zero-based index in the visible runner list.
- `ray_count` or `rays`: number of rays to generate.
- `max_ray_count` or `max_rays_traced`: maximum traced ray count.
- `max_threads` or `threads`: worker thread count.
- `seed_value` or `seed`: random seed.
- `sun_shape` or `include_sun_shape_errors`: enable sun-shape errors.
- `optical_errors` or `include_optical_errors`: enable optical errors.
- `point_focus_system`: enable point-focus-system mode.

Example:

```js
sim.start({
  runner: "embree",
  ray_count: 100000,
  max_ray_count: 1000000,
  seed: 42,
  sun_shape: true,
  optical_errors: true
})
```

## scenes.set_current_index(index)

Select an open scene by zero-based list index. Returns `true` when the scene was selected.

Example:

```js
scenes.set_current_index(0)
```

## scenes.set_current_name(name)

Select the first open scene with the given display name. Returns `true` when a matching scene was selected.

Names are display names and are not guaranteed to be unique.

Example:

```js
scenes.set_current_name("Receiver Study")
```

## scenes.new_blank(name)

Create a new blank scene and select it. Empty names default to `Untitled`. Returns `true` when the request was accepted.

The `db` global is retargeted immediately, so later `db` calls in the same script operate on the new blank scene.

Example:

```js
scenes.new_blank("Generated Scene")
const element = db.create()
```

## scenes.new_from_file(relative_path, name_override)

Load a `.stinput` or `.json` scene file relative to the script working directory. Absolute paths and paths that escape the working directory are rejected. Returns `true` when the load was started.

Loading uses the normal asynchronous file loader, so the loaded scene is not available to later `db` calls in the same script execution.

Example:

```js
scenes.new_from_file("inputs/example.stinput", "Imported Example")
```

## scenes.export_json(relative_path)

Export the current scene as JSON to a path relative to the script working directory. Absolute paths and paths that escape the working directory are rejected. Returns `true` when the export succeeds.

Example:

```js
scenes.export_json("exports/current_scene.json")
```
