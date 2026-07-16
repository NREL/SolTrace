---
vec3

Create a three-component vector.

Inputs: either vec3(value) to expand one scalar to [value, value, value], or vec3(x, y, z). Output: JSON array [x, y, z].

Example: db.vec3(1, 2, 3)
---
vec3_add

Add two vectors component-wise.

Inputs: a and b are vectors as [x, y, z] or objects {x, y, z}; scalar values are also accepted and expanded to all components. Output: JSON array [x, y, z].

Example: db.vec3_add([1, 2, 3], [4, 5, 6])
---
vec3_sub

Subtract one vector from another component-wise.

Inputs: a and b are vectors as [x, y, z] or objects {x, y, z}; scalar values are also accepted and expanded to all components. Output: JSON array [x, y, z].

Example: db.vec3_sub([4, 5, 6], [1, 2, 3])
---
vec3_scale

Scale a vector by a scalar.

Inputs: value is a vector as [x, y, z] or {x, y, z}; scalar values are expanded to all components. scale is a number. Output: JSON array [x, y, z].

Example: db.vec3_scale([1, 2, 3], 10)
---
vec3_dot

Compute the dot product of two vectors.

Inputs: a and b are vectors as [x, y, z] or {x, y, z}; scalar values are expanded to all components. Output: number.

Example: db.vec3_dot([1, 0, 0], [0, 1, 0])
---
vec3_cross

Compute the cross product of two vectors.

Inputs: a and b are vectors as [x, y, z] or {x, y, z}. Output: JSON array [x, y, z]. Invalid inputs return an empty array.

Example: db.vec3_cross([1, 0, 0], [0, 1, 0])
---
vec3_length

Compute vector length.

Input: value is a vector as [x, y, z] or {x, y, z}; scalar values are expanded to all components. Output: number.

Example: db.vec3_length([3, 4, 0])
---
vec3_distance

Compute distance between two vectors.

Inputs: a and b are vectors as [x, y, z] or {x, y, z}. Output: number.

Example: db.vec3_distance([0, 0, 0], [3, 4, 0])
---
vec3_normalize

Normalize a vector.

Input: value is a vector as [x, y, z] or {x, y, z}. Output: JSON array [x, y, z]. Zero-length or invalid inputs return an empty array.

Example: db.vec3_normalize([0, 0, 10])
---
quat

Create a quaternion.

Inputs: w, x, y, z numbers. Output: JSON array [w, x, y, z].

Example: db.quat(1, 0, 0, 0)
---
quat_identity

Create the identity quaternion.

Inputs: none. Output: JSON array [1, 0, 0, 0].

Example: db.quat_identity()
---
quat_from_axis_angle

Create a quaternion from an axis and angle.

Inputs: axis is a vector as [x, y, z] or {x, y, z}; degrees is an angle in degrees. Output: JSON array [w, x, y, z]. Invalid axis returns an empty array.

Example: db.quat_from_axis_angle([0, 0, 1], 45)
---
quat_mul

Multiply two quaternions.

Inputs: a and b are quaternions as [w, x, y, z] or {w, x, y, z}. Output: JSON array [w, x, y, z].

Example: db.quat_mul(db.quat_identity(), db.quat_from_axis_angle([0, 0, 1], 45))
---
quat_conjugate

Conjugate a quaternion.

Input: value is [w, x, y, z] or {w, x, y, z}. Output: JSON array [w, x, y, z].

Example: db.quat_conjugate([0.707, 0, 0, 0.707])
---
quat_inverse

Invert a quaternion.

Input: value is [w, x, y, z] or {w, x, y, z}. Output: JSON array [w, x, y, z]. Zero-length or invalid inputs return an empty array.

Example: db.quat_inverse([0.707, 0, 0, 0.707])
---
quat_normalize

Normalize a quaternion.

Input: value is [w, x, y, z] or {w, x, y, z}. Output: JSON array [w, x, y, z]. Zero-length or invalid inputs return an empty array.

Example: db.quat_normalize([2, 0, 0, 0])
---
quat_rotate_vec3

Rotate a vector by a quaternion.

Inputs: rotation is [w, x, y, z] or {w, x, y, z}; value is [x, y, z] or {x, y, z}. Output: JSON array [x, y, z].

Example: db.quat_rotate_vec3(db.quat_from_axis_angle([0, 0, 1], 90), [1, 0, 0])
---
get_ray_source

Read the current ray source position, shape, and ray generation mode.

Inputs: none. Output: object {type, position, generation, shape, sigma, half_width, csr, user_distribution}. position is [x, y, z]. generation is random or halton. shape is gaussian, pillbox, limb_darkened, buie_csr, user_defined, none, or unknown. user_distribution is an array of {angle, intensity}.

Example: const source = db.get_ray_source()
---
set_ray_source

Update the ray source from an object.

Input: source object may contain position [x, y, z], generation random|halton, shape gaussian|pillbox|limb_darkened|buie_csr|user_defined, sigma, half_width, csr, and user_distribution. user_distribution accepts [{angle, intensity}, ...] or [[angle, intensity], ...]. Output: none.

Example: db.set_ray_source({ position: [0, 0, 1000], generation: "halton", shape: "pillbox", half_width: 4.65 })
---
set_sun_direction

Set the sun direction from a three-component vector.

Input: direction is [x, y, z] or {x, y, z}. The vector is normalized before storage. Output: none.

Example: db.set_sun_direction([0, 0, -1])
---
set_sun_position

Set the sun position from a three-component vector.

Input: position is [x, y, z] or {x, y, z}. Output: none.

Example: db.set_sun_position([0, 0, 1000])
---
set_sun_shape

Update only the sun shape and shape parameters.

Input: shape object may contain shape, sigma, half_width, csr, and user_distribution. Shapes are gaussian, pillbox, limb_darkened, buie_csr, user_defined, none, or custom as an alias for user_defined. Output: none.

Example: db.set_sun_shape({ shape: "gaussian", sigma: 2.8 })
---
get_all_elements

Obtain all scene elements.

Inputs: none. Output: array of Entity handles for all entities with an ElementComponent. Use get_identity, get_transform, get_material_of, and get_geometry_of to inspect each element.

Example: const elements = db.get_all_elements()
---
create

Create a new CSP element. Requires geometry and material to be useful.

Inputs: none. Output: Entity handle for a new scene element with an identity transform. Assign a material and geometry before simulation.

Example: const element = db.create()
---
destroy

Destroy an entity or remove a material, geometry, or tag group.

Input: entity handle returned by this API. For material and geometry groups, related assignments are cleaned up through the database API. For scene elements, children are promoted to top-level elements. Output: none.

Example: db.destroy(element)
---
valid

Check whether an entity handle still points to an existing entity.

Input: entity handle. Output: boolean.

Example: if (db.valid(element)) { db.set_identity(element, "Receiver") }
---
get_identity

Get an entity display name.

Input: entity handle. Output: string. Invalid entities return an empty string.

Example: const name = db.get_identity(element)
---
set_identity

Set an entity display name.

Inputs: entity handle and name string. Output: none.

Example: db.set_identity(element, "Heliostat 1")
---
get_invisible

Check whether an entity is hidden in visualization.

Input: entity handle. Output: boolean. Invalid entities return false.

Example: const hidden = db.get_invisible(element)
---
set_invisible

Set whether an entity is hidden in visualization.

Inputs: entity handle and invisible boolean. Output: none.

Example: db.set_invisible(element, true)
---
get_transform

Read an entity transform.

Input: entity handle. Output: object {position, rotation}, where position is [x, y, z] and rotation is quaternion [w, x, y, z]. Invalid entities return an empty object.

Example: const transform = db.get_transform(element)
---
set_transform

Patch an entity transform.

Inputs: entity handle and transform object. transform may contain position [x, y, z] or {x, y, z}, and rotation [w, x, y, z] or {w, x, y, z}. Output: none.

Example: db.set_transform(element, { position: [1, 2, 3], rotation: db.quat_identity() })
---
get_all_materials

Obtain all available materials, identified by internal IDs.

Inputs: none. Output: array of Entity handles for material groups.

Example: const materials = db.get_all_materials()
---
create_material

Create a new material group.

Inputs: none. Output: Entity handle for the created material group. Use set_identity and set_material_properties to configure it.

Example: const material = db.create_material()
---
get_material_properties

Read material optical properties.

Input: material Entity handle. Output: object {front, back}. Each side contains my_type, error_distribution_type, transmissivity, reflectivity, slope_error, specularity_error, refraction_index_front, and refraction_index_back.

Example: const props = db.get_material_properties(material)
---
set_material_properties

Patch material optical properties.

Inputs: material Entity handle and properties object {front, back}. Each side may contain my_type, error_distribution_type, transmissivity, reflectivity, slope_error, specularity_error, refraction_index_front, and refraction_index_back. Output: none.

Example: db.set_material_properties(material, { front: { reflectivity: 0.95 }, back: { reflectivity: 0.95 } })
---
remove_material

Remove a material group.

Input: material Entity handle. Output: none.

Example: db.remove_material(material)
---
get_all_geometries

Obtain all available geometry groups.

Inputs: none. Output: array of Entity handles for geometry groups.

Example: const geometries = db.get_all_geometries()
---
create_geometry

Create a new geometry group.

Inputs: none. Output: Entity handle for the created geometry group. Use set_identity and set_geometry_properties to configure it.

Example: const geometry = db.create_geometry()
---
get_geometry_properties

Read geometry properties.

Input: geometry Entity handle. Output: object that may contain aperture and surface JSON objects in the same shape used by the simulation data format.

Example: const props = db.get_geometry_properties(geometry)
---
set_geometry_properties

Patch geometry properties.

Inputs: geometry Entity handle and properties object. The object may contain aperture and surface JSON objects in the same shape used by the simulation data format. Output: none.

Example: db.set_geometry_properties(geometry, { aperture: { type: "rectangle", width: 1, height: 1 } })
---
remove_geometry

Remove a geometry group.

Input: geometry Entity handle. Output: none.

Example: db.remove_geometry(geometry)
---
get_material_of

Get the material assigned to a scene element.

Input: scene element Entity handle. Output: material Entity handle, or an invalid handle when unassigned.

Example: const material = db.get_material_of(element)
---
set_material_of

Assign a material to a scene element.

Inputs: scene element Entity handle and material Entity handle. Output: none.

Example: db.set_material_of(element, material)
---
get_geometry_of

Get the geometry assigned to a scene element.

Input: scene element Entity handle. Output: geometry Entity handle, or an invalid handle when unassigned.

Example: const geometry = db.get_geometry_of(element)
---
set_geometry_of

Assign geometry to a scene element.

Inputs: scene element Entity handle and geometry Entity handle. Output: none.

Example: db.set_geometry_of(element, geometry)
---
get_text_content

Read a UTF-8 text file from the script working directory.

Input: relative_path is a path relative to the script working directory. Absolute paths and paths that escape the working directory are rejected. Output: string content, or an empty string if the file cannot be opened.

Example: const csv = db.get_text_content("inputs/heliostats.csv")
---
get_json_content

Read a JSON file from the script working directory.

Input: relative_path is a path relative to the script working directory. Absolute paths and paths that escape the working directory are rejected. Output: parsed JSON object or array, or an empty value if the file cannot be opened or parsed.

Example: const layout = db.get_json_content("inputs/layout.json")
