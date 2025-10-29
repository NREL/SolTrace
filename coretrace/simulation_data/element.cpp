
#include "element.hpp"

#include <math.h>

#include "constants.hpp"

namespace SolTrace::Data {

// ElementContainer ElementBase::empty_container;

ElementBase::ElementBase() : Element(),
                             coordinates_initialized(true),
                             active(true),
                             virtual_flag(false),
                             my_id(ELEMENT_ID_UNASSIGNED),
                             stage(-1),
                             zrot(0.0),
                             reference_element(nullptr)
{
    // Default local coordinates to match with the reference coordinates
    this->aim.set_values(0.0, 0.0, 1.0);
    this->origin.zero();

    this->euler_angles.zero();
    this->reference_to_local.identity();
    this->local_to_reference.identity();

    return;
}

ElementBase::~ElementBase()
{
    this->reference_element = nullptr;
    return;
}

// ElementContainer::iterator ElementBase::get_iterator()
// {
//     return empty_container.get_iterator();
// }

// ElementContainer::const_iterator ElementBase::get_const_iterator()
// {
//     return empty_container.get_const_iterator();
// }

Vector3d ElementBase::get_origin_stage() const
{
    Vector3d origin_stage;
    auto ref_el = this->reference_element;
    if (this->is_stage())
    {
        // This is the stage element so the stage origin is zero.
        origin_stage.zero();
    }
    else if (ref_el == nullptr)
    {
        // We hit the end of the chain without finding a stage element.
        // Assume we are not using stages and return the global origin.
        origin_stage = this->origin;
    }
    else
    {
        ref_el->convert_local_to_stage(origin_stage, this->origin);
    }
    return origin_stage;
}

Vector3d ElementBase::get_origin_global() const
{
    Vector3d origin_global;
    auto ref_el = this->reference_element;
    if (ref_el == nullptr)
    {
        origin_global = this->origin;
    }
    else
    {
        ref_el->convert_local_to_global(origin_global, this->origin);
    }
    return origin_global;
}

Vector3d ElementBase::get_aim_vector_stage() const
{
    Vector3d aim_stage;
    if (this->reference_element == nullptr)
    {
        aim_stage = this->aim;
    }
    else
    {
        this->reference_element->convert_local_to_stage(
            aim_stage, this->aim);
    }
    return aim_stage;
}

Vector3d ElementBase::get_aim_vector_global() const
{
    Vector3d aim_global;
    if (this->reference_element == nullptr)
    {
        aim_global = this->aim;
    }
    else
    {
        this->reference_element->convert_local_to_global(
            aim_global, this->aim);
    }
    return aim_global;
}

Matrix3d ElementBase::get_reference_to_local() const
{
    return this->reference_to_local;
}

Matrix3d ElementBase::get_stage_to_local() const
{
    Matrix3d stage_to_local;
    if (this->is_stage())
    {
        // stage_to_local.set_value(0, 0, 1.0);
        // stage_to_local.set_value(1, 1, 1.0);
        // stage_to_local.set_value(2, 2, 1.0);
        stage_to_local.identity();
    }
    else if (this->reference_element == nullptr)
    {
        stage_to_local = this->reference_to_local;
    }
    else
    {
        Matrix3d R = this->reference_element->get_stage_to_local();
        matrix_matrix_product(this->reference_to_local, R, stage_to_local);
    }
    return stage_to_local;
}

Matrix3d ElementBase::get_global_to_local() const
{
    Matrix3d global_to_local;
    if (this->reference_element == nullptr)
    {
        global_to_local = this->reference_to_local;
    }
    else
    {
        Matrix3d R = this->reference_element->get_global_to_local();
        matrix_matrix_product(this->reference_to_local, R, global_to_local);
    }
    return global_to_local;
}

Matrix3d ElementBase::get_local_to_reference() const
{
    return this->local_to_reference;
}

Matrix3d ElementBase::get_local_to_stage() const
{
    Matrix3d local_to_stage;
    if (this->is_stage())
    {
        local_to_stage.set_value(0, 0, 1.0);
        local_to_stage.set_value(1, 1, 1.0);
        local_to_stage.set_value(2, 2, 1.0);
    }
    else if (this->reference_element == nullptr)
    {
        local_to_stage = this->local_to_reference;
    }
    else
    {
        Matrix3d R = this->reference_element->get_local_to_stage();
        matrix_matrix_product(R, this->local_to_reference, local_to_stage);
    }
    return local_to_stage;
}

Matrix3d ElementBase::get_local_to_global() const
{
    Matrix3d local_to_global;
    if (this->reference_element == nullptr)
    {
        local_to_global = this->local_to_reference;
    }
    else
    {
        Matrix3d R = this->reference_element->get_local_to_global();
        matrix_matrix_product(R, this->local_to_reference, local_to_global);
    }
    return local_to_global;
}

int ElementBase::compute_coordinate_rotations()
{
    int sts = 0;

    if (!this->coordinates_initialized)
    {
        Vector3d dr;
        vector_add(1.0, this->aim, -1.0, this->origin, dr);
        make_unit_vector(dr);

        this->euler_angles[0] = atan2(dr[0], dr[2]);
        this->euler_angles[1] = asin(dr[1]);
        this->euler_angles[2] = this->zrot * D2R;

        compute_transform_matrices(this->euler_angles,
                                   this->reference_to_local,
                                   this->local_to_reference);

        this->coordinates_initialized = true;
    }

    return sts;
}

int ElementBase::set_reference_frame_geometry(const Vector3d &origin,
                                              const Vector3d &aim,
                                              double zrot)
{
    this->coordinates_initialized = false;
    this->origin = origin;
    this->aim = aim;
    this->zrot = zrot;
    return this->compute_coordinate_rotations();
}

int ElementBase::convert_reference_to_local(Vector3d &local,
                                            const Vector3d &ref)
{
    Vector3d temp;
    vector_add(1.0, ref, -1.0, this->origin, temp);
    matrix_vector_product(this->reference_to_local, temp, local);
    return 0;
}

int ElementBase::convert_stage_to_local(Vector3d &local,
                                        const Vector3d &stage)
{
    if (this->is_stage())
    {
        // We are in the stage coordinate frame so the local coordinates
        // are the stage coordinates
        local = stage;
    }
    else if (this->reference_element == nullptr)
    {
        // No stage frame found so assume we are not using stages
        // and return the same answer as convert_global_to_local
        this->convert_reference_to_local(local, stage);
    }
    else
    {
        Vector3d ref;
        this->reference_element->convert_stage_to_local(ref, stage);
        convert_reference_to_local(local, ref);
    }
    return 0;
}

int ElementBase::convert_global_to_local(Vector3d &local,
                                         const Vector3d &global)
{
    auto ref_el = this->reference_element;
    if (ref_el == nullptr)
    {
        // We are at the global coordinate frame
        this->convert_reference_to_local(local, global);
    }
    else
    {
        Vector3d ref;
        ref_el->convert_global_to_local(ref, global);
        this->convert_reference_to_local(local, ref);
    }
    return 0;
}

int ElementBase::convert_local_to_reference(Vector3d &ref,
                                            const Vector3d &local)
{
    Vector3d temp;
    matrix_vector_product(this->local_to_reference, local, temp);
    vector_add(1.0, temp, 1.0, this->origin, ref);
    return 0;
}

int ElementBase::convert_local_to_stage(Vector3d &stage,
                                        const Vector3d &local)
{
    if (this->is_stage())
    {
        // We are in the stage coordinate frame so the local coordinates
        // are the stage coordinates
        stage = local;
    }
    else if (this->reference_element == nullptr)
    {
        // No stage has been found so assume we are not using stages
        // and convert to global coordinates
        this->convert_local_to_reference(stage, local);
        // TODO: This should probably return something other than 0...
    }
    else
    {
        Vector3d ref;
        this->convert_local_to_reference(ref, local);
        this->reference_element->convert_local_to_stage(stage, ref);
    }
    return 0;
}

int ElementBase::convert_local_to_global(Vector3d &global,
                                         const Vector3d &local)
{
    auto ref_el = this->reference_element;
    if (ref_el == nullptr)
    {
        this->convert_local_to_reference(global, local);
    }
    else
    {
        Vector3d ref;
        this->convert_local_to_reference(ref, local);
        ref_el->convert_local_to_global(global, ref);
    }
    return 0;
}

int ElementBase::convert_global_to_reference(Vector3d &ref,
                                             const Vector3d &global)
{
    if (this->reference_element == nullptr)
    {
        ref = global;
        return 0;
    }
    else
    {
        return this->reference_element->convert_global_to_local(ref, global);
    }
}

int ElementBase::convert_reference_to_global(Vector3d &global,
                                             const Vector3d &ref)
{
    if (this->reference_element == nullptr)
    {
        global = ref;
        return 0;
    }
    else
    {
        return this->reference_element->convert_local_to_global(global, ref);
    }
}

int ElementBase::convert_vector_reference_to_local(Vector3d &local,
                                                   const Vector3d &ref)
{
    matrix_vector_product(this->reference_to_local, ref, local);
    return 0;
}

int ElementBase::convert_vector_stage_to_local(Vector3d &local,
                                               const Vector3d &stage)
{
    if (this->is_stage())
    {
        // We are in the stage coordinate frame so the local coordinates
        // are the stage coordinates
        local = stage;
    }
    else if (this->reference_element == nullptr)
    {
        // No stage frame found so assume we are not using stages
        // and return the same answer as convert_global_to_local
        this->convert_vector_reference_to_local(local, stage);
    }
    else
    {
        Vector3d ref;
        this->reference_element->convert_vector_stage_to_local(ref, stage);
        convert_vector_reference_to_local(local, ref);
    }
    return 0;
}

int ElementBase::convert_vector_global_to_local(Vector3d &local,
                                                const Vector3d &global)
{
    auto ref_el = this->reference_element;
    if (ref_el == nullptr)
    {
        // We are at the global coordinate frame
        this->convert_vector_reference_to_local(local, global);
    }
    else
    {
        Vector3d ref;
        ref_el->convert_vector_global_to_local(ref, global);
        this->convert_vector_reference_to_local(local, ref);
    }
    return 0;
}

int ElementBase::convert_vector_local_to_reference(Vector3d &ref,
                                                   const Vector3d &local)
{
    matrix_vector_product(this->local_to_reference, local, ref);
    return 0;
}

int ElementBase::convert_vector_local_to_stage(Vector3d &stage,
                                               const Vector3d &local)
{
    if (this->is_stage())
    {
        // We are in the stage coordinate frame so the local coordinates
        // are the stage coordinates
        stage = local;
    }
    else if (this->reference_element == nullptr)
    {
        // No stage has been found so assume we are not using stages
        // and convert to global coordinates
        this->convert_vector_local_to_reference(stage, local);
        // TODO: This should probably return something other than 0...
    }
    else
    {
        Vector3d ref;
        this->convert_vector_local_to_reference(ref, local);
        this->reference_element->convert_vector_local_to_stage(stage, ref);
    }
    return 0;
}

int ElementBase::convert_vector_local_to_global(Vector3d &global,
                                                const Vector3d &local)
{
    auto ref_el = this->reference_element;
    if (ref_el == nullptr)
    {
        this->convert_vector_local_to_reference(global, local);
    }
    else
    {
        Vector3d ref;
        this->convert_vector_local_to_reference(ref, local);
        ref_el->convert_vector_local_to_global(global, ref);
    }
    return 0;
}

int ElementBase::convert_vector_global_to_reference(Vector3d &ref,
                                                    const Vector3d &global)
{
    if (this->reference_element == nullptr)
    {
        ref = global;
        return 0;
    }
    else
    {
        return this->reference_element->convert_vector_global_to_local(
            ref, global);
    }
}

int ElementBase::convert_vector_reference_to_global(Vector3d &global,
                                                    const Vector3d &ref)
{
    if (this->reference_element == nullptr)
    {
        global = ref;
        return 0;
    }
    else
    {
        return this->reference_element->convert_vector_local_to_global(
            global, ref);
    }
}

// const OpticalProperties & ElementBase::get_optical_properties() const
// {
//     return this->optics;
// }

// void ElementBase::set_optical_properties(const OpticalProperties &op)
// {
//     this->optics = op;
// }

void ElementBase::enforce_user_fields_set() const
{
    return;
}

} // namespace SolTrace::Data
