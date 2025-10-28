/**
 * @file element.hpp
 * @brief Base element class and element management
 *
 * Defines the base ElementBase class and element management utilities.
 * All optical elements in SolTrace derive from this base class, providing
 * common functionality for positioning, orientation, and optical properties.
 *
 * @defgroup elements Optical Elements
 * @{
 */

#ifndef SOLTRACE_ELEMENT_H
#define SOLTRACE_ELEMENT_H

// #include <memory>
#include <string>

#include "aperture.hpp"
#include "constants.hpp"
#include "container.hpp"
#include "optical_properties.hpp"
#include "ray_source.hpp"
#include "surface.hpp"
#include "vector3d.hpp"

namespace SolTrace::Data {

using element_id = std::int_fast64_t;
const element_id ELEMENT_ERROR = -1;
const element_id ELEMENT_ID_UNASSIGNED = -2;
const element_id ELEMENT_ALREADY_REGISTERED = -3;
const element_id ELEMENT_INVALID_SETUP = -4;
const element_id ELEMENT_NULL = -5;

// Forward declaration of the Element class so we can define ElementContainer
class Element;

using ElementContainer = Container<element_id, Element>;
using element_ptr = ElementContainer::value_pointer;

class Element
{
public:
  static bool is_success(element_id id)
  {
    return id >= 0;
  }

  Element() {};
  virtual ~Element() {};

  // Accessors for any element
  /// @brief Disable the element for ray tracing
  virtual void disable() const = 0;
  /// @brief Enable the element for ray tracing
  virtual void enable() const = 0;
  /// @brief Check whether element is enabled/disabled
  /// @return true if enabled, false otherwise
  virtual bool is_enabled() const = 0;

  /// @brief Check whether the element is a CompositeElement
  /// @return true if CompositeElement, false otherwise
  virtual bool is_composite() const = 0;
  /// @brief Check whether the element is a SingleElement
  /// @return true if SingleElement, false otherwise
  virtual bool is_single() const = 0;
  /// @brief Check whether this element is a StageElement
  /// @return true if StageElement, false otherwise
  virtual bool is_stage() const = 0;
  /// @brief Check whether the element is a VirtualElement
  /// @return true if VirtualElement, false otherwise
  virtual bool is_virtual() const = 0;
  virtual void mark_virtual() const = 0;
  virtual void unmark_virtual() const = 0;

  /// @brief Get the element id assigned when registered with SimulationData
  /// @return id if registered with SimulationData, ELEMENT_ID_UNASSIGNED if not
  virtual element_id get_id() const = 0;

  /**
   * @brief Get the stage number this element belongs to
   * @return Stage number
   */
  virtual int_fast64_t get_stage() = 0;

  /**
   * @brief Get the element name
   * @return Reference to element name string
   */
  virtual const std::string &get_name() const = 0;

  /**
   * @brief Set the element name
   * @param name New name for the element
   */
  virtual void set_name(const std::string &name) = 0;

  /****************************************************************************
   * NOTE: For all coordinate functions below, the term "reference coordinates"
   * means the coordinate frame immediately above this elements. If Element
   * is a subelement of a CompositeElement the reference coordinates are the
   * CompositeElements even if that CompositeElement is then stored in a stage.
   ***************************************************************************/

  /**
   * @brief Get origin position in reference coordinates
   * @return Origin position vector in reference frame
   */
  virtual Vector3d get_origin_ref() const = 0;

  /**
   * @brief Get origin position in stage coordinates
   * @return Origin position vector in stage frame
   */
  virtual Vector3d get_origin_stage() const = 0;

  /**
   * @brief Get origin position in global coordinates
   * @return Origin position vector in global frame
   */
  virtual Vector3d get_origin_global() const = 0;

  /**
   * @brief Set origin position (always relative to reference coordinates)
   * @param origin New origin position vector
   */
  virtual void set_origin(const Vector3d &) = 0;

  /**
   * @brief Set origin position (always relative to reference coordinates)
   * @param x X coordinate
   * @param y Y coordinate
   * @param z Z coordinate
   */
  virtual void set_origin(double, double, double) = 0;

  // virtual const Vector3d &get_global_origin() const = 0;
  // virtual void set_global_origin(const Vector3d &) = 0;

  /**
   * @brief Get aim vector in reference coordinates
   * @return Aim direction vector in reference frame
   */
  virtual Vector3d get_aim_vector_ref() const = 0;
  virtual Vector3d get_aim_vector_stage() const = 0;
  virtual Vector3d get_aim_vector_global() const = 0;
  // Always the aim vector with respect the reference coordinates
  virtual void set_aim_vector(const Vector3d &) = 0;
  virtual void set_aim_vector(double, double, double) = 0;
  // Always the Euler angles with respect the reference coordinates
  virtual const Vector3d &get_euler_angles() const = 0;
  // Always the ZRot with respect to the reference coordinates
  virtual double get_zrot() const = 0;
  virtual void set_zrot(double) = 0;
  virtual double get_zrot_radians() const = 0;
  virtual void set_zrot_radians(double) = 0;

  virtual Matrix3d get_reference_to_local() const = 0;
  virtual Matrix3d get_stage_to_local() const = 0;
  virtual Matrix3d get_global_to_local() const = 0;
  virtual Matrix3d get_local_to_reference() const = 0;
  virtual Matrix3d get_local_to_stage() const = 0;
  virtual Matrix3d get_local_to_global() const = 0;

  // virtual const Vector3d &get_upper_bounding_box() const = 0;
  // virtual const Vector3d &get_lower_bounding_box() const = 0;

  // Accessors for SingleElements
  virtual const aperture_ptr get_aperture() const = 0;
  virtual aperture_ptr get_aperture() = 0;
  virtual void set_aperture(aperture_ptr) = 0;
  virtual const surface_ptr get_surface() const = 0;
  virtual surface_ptr get_surface() = 0;
  virtual void set_surface(surface_ptr) = 0;

  // virtual const OpticalProperties &get_optical_properties() const = 0;
  // virtual void set_optical_properties(const OpticalProperties &) = 0;

  virtual const OpticalProperties *get_front_optical_properties() const = 0;
  virtual OpticalProperties *get_front_optical_properties() = 0;
  virtual void set_front_optical_properties(const OpticalProperties &) = 0;

  virtual const OpticalProperties *get_back_optical_properties() const = 0;
  virtual OpticalProperties *get_back_optical_properties() = 0;
  virtual void set_back_optical_properties(const OpticalProperties &) = 0;

  // Accessors for CompositeElements
  virtual uint_fast64_t get_number_of_elements() const = 0;
  // virtual ElementContainer::iterator get_iterator() = 0;
  // virtual ElementContainer::const_iterator get_const_iterator() = 0;
  // virtual bool is_at_end(ElementContainer::iterator iter) = 0;
  // virtual bool is_at_end(ElementContainer::const_iterator iter) = 0;

  // Coordinate transformation routines
  virtual int set_reference_frame_geometry(const Vector3d &origin,
                                           const Vector3d &aim,
                                           double zrot) = 0;

  /****************************************************************
   * These are point coordinate conversion routines. They convert a point
   * from one coordinate system to another. Assumes that the element
   * hierarchy is set and that `compute_coordinate_rotations` has been
   * called.
   ****************************************************************/

  // Convert `ref` to local coordinates and store the result in `local`
  virtual int convert_reference_to_local(Vector3d &local,
                                         const Vector3d &ref) = 0;
  // Convert `stage` to local coordinates and store the result in `local`
  virtual int convert_stage_to_local(Vector3d &local,
                                     const Vector3d &stage) = 0;
  // Convert `global` to local coordinates and store the result in `local`
  virtual int convert_global_to_local(Vector3d &local,
                                      const Vector3d &global) = 0;
  // Convert `local` to reference coordinates and store the result in `ref`
  virtual int convert_local_to_reference(Vector3d &ref,
                                         const Vector3d &local) = 0;
  // Convert `local` to stage coordinates and store the result in `stage`
  virtual int convert_local_to_stage(Vector3d &stage,
                                     const Vector3d &local) = 0;
  // Convert `local` to global coordinates and store the result in `global`
  virtual int convert_local_to_global(Vector3d &global,
                                      const Vector3d &local) = 0;

  // Convert global coordinates to reference coordinates
  virtual int convert_global_to_reference(Vector3d &ref,
                                          const Vector3d &global) = 0;
  // Convert reference coordinates to global
  virtual int convert_reference_to_global(Vector3d &global,
                                          const Vector3d &ref) = 0;

  /****************************************************************
   * These are vector coordinate conversion routines. They convert a
   * vector (i.e. origin is always the same) from one coordinate
   * system to another. Assumes that the element hierarchy is set
   * and that `compute_coordinate_rotations` has been called.
   ****************************************************************/

  // Convert `ref` to local coordinates and store the result in `local`
  virtual int convert_vector_reference_to_local(Vector3d &local,
                                                const Vector3d &ref) = 0;
  // Convert `stage` to local coordinates and store the result in `local`
  virtual int convert_vector_stage_to_local(Vector3d &local,
                                            const Vector3d &stage) = 0;
  // Convert `global` to local coordinates and store the result in `local`
  virtual int convert_vector_global_to_local(Vector3d &local,
                                             const Vector3d &global) = 0;
  // Convert `local` to reference coordinates and store the result in `ref`
  virtual int convert_vector_local_to_reference(Vector3d &ref,
                                                const Vector3d &local) = 0;
  // Convert `local` to stage coordinates and store the result in `stage`
  virtual int convert_vector_local_to_stage(Vector3d &stage,
                                            const Vector3d &local) = 0;
  // Convert `local` to global coordinates and store the result in `global`
  virtual int convert_vector_local_to_global(Vector3d &global,
                                             const Vector3d &local) = 0;

  // Convert global coordinates to reference coordinates
  virtual int convert_vector_global_to_reference(Vector3d &ref,
                                                 const Vector3d &global) = 0;
  // Convert reference coordinates to global
  virtual int convert_vector_reference_to_global(Vector3d &global,
                                                 const Vector3d &ref) = 0;

  // Other routines
  // Computes necessary coordinate transformation data. Expects
  // that the element hierarchy is set above this element.
  virtual int compute_coordinate_rotations() = 0;

  // WARNING: The below Accessors should be used with EXTREME caution!!!
  // These are used by other classes to set things up correctly and
  // not meant for the casual user. You have been warned!

  /// @brief Set the element id--used by SimulationData
  /// @param id id assigned and to set
  virtual void set_id(element_id id) = 0;
  virtual void set_reference_element(Element *reference) = 0;
  virtual void set_stage(int_fast64_t stage) = 0;

  // WARNING: The below Accessors should be used with care. They set
  // values that are set automatically -- these are here just in case...
  virtual void set_euler_angles(const Vector3d &) = 0;
  virtual void set_euler_angles(double, double, double) = 0;
  virtual void set_reference_to_local(const Matrix3d &) = 0;
  virtual void set_local_to_reference(const Matrix3d &) = 0;

  // Check that all required fields have been set
  virtual void enforce_user_fields_set() const = 0;

protected:
  // virtual int set_bounding_box() = 0;

private:
};

class ElementBase : public Element
{
public:
  ElementBase();
  // ElementBase(const Vector3d &origin, const Vector3d &aim);
  virtual ~ElementBase();

  virtual inline void disable() const override { this->active = false; }
  virtual inline void enable() const override { this->active = true; }
  virtual bool is_enabled() const override { return this->active; }

  virtual bool is_composite() const override { return false; }
  virtual bool is_single() const override { return false; }
  virtual bool is_stage() const override { return false; }

  virtual bool is_virtual() const override { return this->virtual_flag; }
  virtual void mark_virtual() const override {this->virtual_flag = true;}
  virtual void unmark_virtual() const override {this->virtual_flag = false;}

  // virtual ElementContainer::iterator get_iterator();
  // virtual ElementContainer::const_iterator get_const_iterator();
  // virtual bool is_at_end(ElementContainer::iterator iter) { return true; }
  // virtual bool is_at_end(ElementContainer::const_iterator iter) { return true; }
  virtual void set_reference_element(Element *reference) override
  {
    this->reference_element = reference;
  }

  virtual uint_fast64_t get_number_of_elements() const override { return 1; }

  virtual element_id get_id() const override { return this->my_id; }
  virtual void set_id(element_id id) override
  {
    this->my_id = id;
    return;
  }

  virtual int_fast64_t get_stage() override { return this->stage; }
  virtual void set_stage(int_fast64_t stage) override { this->stage = stage; }

  virtual const std::string &get_name() const override
  {
    return this->my_name;
  }
  virtual void set_name(const std::string &name) override
  {
    this->my_name = name;
  }

  virtual Vector3d get_origin_ref() const override { return this->origin; }
  virtual Vector3d get_origin_stage() const override;
  virtual Vector3d get_origin_global() const override;
  virtual void set_origin(const Vector3d &point) override
  {
    this->coordinates_initialized = false;
    this->origin = point;
    return;
  }
  virtual void set_origin(double x, double y, double z) override
  {
    this->coordinates_initialized = false;
    this->origin.set_values(x, y, z);
    return;
  }
  virtual Vector3d get_aim_vector_ref() const override { return this->aim; }
  virtual Vector3d get_aim_vector_stage() const override;
  virtual Vector3d get_aim_vector_global() const override;
  virtual void set_aim_vector(const Vector3d &direction) override
  {
    this->coordinates_initialized = false;
    this->aim = direction;
    return;
  }
  virtual void set_aim_vector(double x, double y, double z) override
  {
    this->coordinates_initialized = false;
    this->aim.set_values(x, y, z);
    return;
  }
  virtual const Vector3d &get_euler_angles() const override
  {
    return this->euler_angles;
  }
  // virtual void set_euler_angles(const Vector3d &angles)
  // {
  //   this->euler_angles = angles;
  //   return;
  // }
  virtual double get_zrot() const override { return this->zrot; }
  virtual void set_zrot(double rot) override
  {
    this->coordinates_initialized = false;
    this->zrot = rot;
    return;
  }

  virtual double get_zrot_radians() const override
  {
    return this->zrot * PI / 180.0;
  }
  virtual void set_zrot_radians(double zrad) override
  {
    this->coordinates_initialized = false;
    this->zrot = zrad * 180.0 / PI;
    return;
  }

  virtual Matrix3d get_reference_to_local() const override;
  virtual Matrix3d get_stage_to_local() const override;
  virtual Matrix3d get_global_to_local() const override;
  virtual Matrix3d get_local_to_reference() const override;
  virtual Matrix3d get_local_to_stage() const override;
  virtual Matrix3d get_local_to_global() const override;

  virtual int compute_coordinate_rotations() override;
  virtual int set_reference_frame_geometry(const Vector3d &origin,
                                           const Vector3d &aim,
                                           double zrot) override;

  // Convert `ref` to local coordinates and store the result in `local`
  virtual int convert_reference_to_local(Vector3d &local,
                                         const Vector3d &ref) override;
  // Convert `stage` to local coordinates and store the result in `local`
  virtual int convert_stage_to_local(Vector3d &local,
                                     const Vector3d &stage) override;
  // Convert `global` to local coordinates and store the result in `local`
  virtual int convert_global_to_local(Vector3d &local,
                                      const Vector3d &global) override;
  // Convert `local` to reference coordinates and store the result in `ref`
  virtual int convert_local_to_reference(Vector3d &ref,
                                         const Vector3d &local) override;
  // Convert `local` to stage coordinates and store the result in `stage`
  virtual int convert_local_to_stage(Vector3d &stage,
                                     const Vector3d &local) override;
  // Convert `local` to global coordinates and store the result in `global`
  virtual int convert_local_to_global(Vector3d &global,
                                      const Vector3d &local) override;

  // Convert global coordinates to reference coordinates
  virtual int convert_global_to_reference(Vector3d &ref,
                                          const Vector3d &global) override;
  // Convert reference coordinates to global
  virtual int convert_reference_to_global(Vector3d &global,
                                          const Vector3d &ref) override;

  // Convert `ref` to local coordinates and store the result in `local`
  virtual int convert_vector_reference_to_local(Vector3d &local,
                                                const Vector3d &ref) override;
  // Convert `stage` to local coordinates and store the result in `local`
  virtual int convert_vector_stage_to_local(Vector3d &local,
                                            const Vector3d &stage) override;
  // Convert `global` to local coordinates and store the result in `local`
  virtual int convert_vector_global_to_local(Vector3d &local,
                                             const Vector3d &global) override;
  // Convert `local` to reference coordinates and store the result in `ref`
  virtual int convert_vector_local_to_reference(Vector3d &ref,
                                                const Vector3d &local) override;
  // Convert `local` to stage coordinates and store the result in `stage`
  virtual int convert_vector_local_to_stage(Vector3d &stage,
                                            const Vector3d &local) override;
  // Convert `local` to global coordinates and store the result in `global`
  virtual int convert_vector_local_to_global(Vector3d &global,
                                             const Vector3d &local) override;

  // Convert global coordinates to reference coordinates
  virtual int convert_vector_global_to_reference(Vector3d &ref,
                                                 const Vector3d &global) override;
  // Convert reference coordinates to global
  virtual int convert_vector_reference_to_global(Vector3d &global,
                                                 const Vector3d &ref) override;

  // WARNING: The below Accessors should be used with care. They set
  // values that are set automatically -- these are here just in case...
  virtual void set_euler_angles(const Vector3d &ea) override
  {
    this->euler_angles = ea;
  }
  virtual void set_euler_angles(double alpha, double beta, double gamma) override
  {
    this->euler_angles.set_values(alpha, beta, gamma);
  }
  virtual void set_reference_to_local(const Matrix3d &rtol) override
  {
    this->reference_to_local = rtol;
  }
  virtual void set_local_to_reference(const Matrix3d &ltor) override
  {
    this->local_to_reference = ltor;
  }

  virtual void enforce_user_fields_set() const override;

protected:
  // TODO: Do these need to be mutable?
  mutable bool active;
  mutable bool virtual_flag;
  mutable element_id my_id;

  bool coordinates_initialized;

  int_fast64_t stage;
  std::string my_name;

  // Location of the origin in the reference coordinate system
  Vector3d origin;
  // Aim vector of element--aligns with the local (positive) z-axis
  Vector3d aim;
  // Rotation about the aim vector to set local x and y axes in degrees (ugh!)
  double zrot;

  Vector3d euler_angles;

  Matrix3d reference_to_local;
  Matrix3d local_to_reference;

  // element_ptr reference_element;
  Element *reference_element;

private:
  // static ElementContainer empty_container;
};

// using element_ptr = typename std::shared_ptr<Element>;
// using ElementContainer = std::map<element_id, element_ptr>;

template <typename C, typename... Args>
inline auto make_element(Args &&...args)
{
  return ElementContainer::make_pointer<C>(std::forward<Args>(args)...);
}

} // namespace SolTrace::Data

/**
 * @}
 */

#endif
