/*
TITLE Flat Mirror Arc With Cylindrical Absorber
DESC Creates a vertical cylindrical absorber and flat mirror facets arranged on a 180 degree ground arc.
DESC The mirror planes stand on the ground and face inward toward the absorber.
PROPERTY mirror_count integer 9 1..=200
PROPERTY arc_radius real 35 1..1000
PROPERTY mirror_width real 6 0.01..
PROPERTY mirror_height real 4 0.01..
PROPERTY absorber_height real 12 0.01..
PROPERTY absorber_radius real 2 0.01..
*/
const absorber_material = db.create_material()
db.set_identity(absorber_material, "Absorber material")
db.set_material_properties(absorber_material, {
    front: {
        my_type: "REFLECTION",
        error_distribution_type: "GAUSSIAN",
        transmissivity: 0.0,
        reflectivity: 0.0,
        slope_error: 0.0,
        specularity_error: 0.0,
        refraction_index_front: 1.0,
        refraction_index_back: 1.0,
    },
    back: {
        my_type: "REFLECTION",
        error_distribution_type: "GAUSSIAN",
        transmissivity: 0.0,
        reflectivity: 0.0,
        slope_error: 0.0,
        specularity_error: 0.0,
        refraction_index_front: 1.0,
        refraction_index_back: 1.0,
    },
})

const mirror_material = db.create_material()
db.set_identity(mirror_material, "Ideal mirror material")
db.set_material_properties(mirror_material, {
    front: {
        my_type: "REFLECTION",
        error_distribution_type: "GAUSSIAN",
        transmissivity: 0.0,
        reflectivity: 1.0,
        slope_error: 0.0,
        specularity_error: 0.0,
        refraction_index_front: 1.0,
        refraction_index_back: 1.0,
    },
    back: {
        my_type: "REFLECTION",
        error_distribution_type: "GAUSSIAN",
        transmissivity: 0.0,
        reflectivity: 0.0,
        slope_error: 0.0,
        specularity_error: 0.0,
        refraction_index_front: 1.0,
        refraction_index_back: 1.0,
    },
})

console.log("Created materials...")

const absorber_geometry = db.create_geometry()
db.set_identity(absorber_geometry, "Cylindrical absorber geometry")
db.set_geometry_properties(absorber_geometry, {
    aperture: {
        aperture_type: "RECTANGLE",
        x_length: absorber_radius * 2.0,
        y_length: absorber_height,
        x_coord: -absorber_radius,
        y_coord: 0.0,
    },
    surface: {
        surface_type: "CYLINDER",
        radius: absorber_radius,
    },
})

const mirror_geometry = db.create_geometry()
db.set_identity(mirror_geometry, "Flat mirror geometry")
db.set_geometry_properties(mirror_geometry, {
    aperture: {
        aperture_type: "RECTANGLE",
        x_length: mirror_width,
        y_length: mirror_height,
        x_coord: -mirror_width / 2.0,
        y_coord: 0.0,
    },
    surface: {
        surface_type: "FLAT",
    },
})

console.log("Created geometry...")

const absorber = db.create()
db.set_identity(absorber, "Cylindrical absorber")
db.set_transform(absorber, {
    position: [0.0, absorber_radius, 0.0],
    rotation: db.quat_from_axis_angle([1.0, 0.0, 0.0], 90.0),
})
db.set_material_of(absorber, absorber_material)
db.set_geometry_of(absorber, absorber_geometry)

const mirrors = []
const denominator = Math.max(1, mirror_count - 1)

for (let i = 0; i < mirror_count; ++i) {
    const fraction = i / denominator
    const angle_degrees = -90.0 + 180.0 * fraction
    const angle_radians = angle_degrees * Math.PI / 180.0
    const x = arc_radius * Math.sin(angle_radians)
    const y = -arc_radius * Math.cos(angle_radians)
    const rotation = db.quat_mul(
        db.quat_from_axis_angle([0.0, 0.0, 1.0], angle_degrees + 180.0),
        db.quat_from_axis_angle([1.0, 0.0, 0.0], 90.0)
    )

    const mirror = db.create()
    db.set_identity(mirror, "Mirror " + String(i + 1))
    db.set_transform(mirror, {
        position: [x, y, 0.0],
        rotation: rotation,
    })
    db.set_material_of(mirror, mirror_material)
    db.set_geometry_of(mirror, mirror_geometry)

    mirrors.push(mirror)
}

console.log("Done!")
