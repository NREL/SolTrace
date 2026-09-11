#pragma once

#include <cstdint>
#include <string>
#include <memory>

#include "vec3d.h"
#include "soltrace_constants.h"
#include "soltrace_type.h"
#include "Surface.h"
#include "Aperture.h"
#include "utils/math_util.h"
#include "shaders/GeometryDataST.h"
#include "shaders/MaterialDataST.h"

namespace SolTrace::Data
{
    class OpticalPropertySet;
}

namespace OptixCSP
{

    class Aperture;
    class CspElementBase
    {
    public:
        CspElementBase();
        virtual ~CspElementBase() = default;

        // Positioning and orientation.
        virtual const Vec3d &get_origin() const = 0;
        virtual void set_origin(const OptixCSP::Vec3d &) = 0;
      virtual const Vec3d &get_aim_point() const = 0;
      virtual void set_aim_point(const Vec3d &o) = 0;

        // virtual const Vec3d& get_euler_angles() const = 0;
        // virtual void set_euler_angles(const Vec3d&) = 0;

        // Bounding box accessors.
        // virtual const Vec3d& get_upper_bounding_box() const = 0;
        // virtual const Vec3d& get_lower_bounding_box() const = 0;

        virtual GeometryDataST toDeviceGeometryData() const = 0;
        virtual MaterialData toDeviceMaterialDataFront() const = 0;
        virtual MaterialData toDeviceMaterialDataBack() const = 0;

        virtual void set_id(const int32_t id) = 0;

    protected:
        // Derived classes must implement bounding box computation.
        // virtual int set_bounding_box() = 0;
    };

    // A concrete implementation of Element that stores data in member variables.
    class CspElement : public CspElementBase
    {
    public:
        CspElement();
        ~CspElement() = default;

        // set and get origin
        const Vec3d &get_origin() const override;
        void set_origin(const Vec3d &o) override;

      const Vec3d &get_aim_point() const override;
      void set_aim_point(const Vec3d &o) override;

        std::shared_ptr<Aperture> get_aperture() const;
        std::shared_ptr<Surface> get_surface() const;
        ApertureType get_aperture_type() const;
        SurfaceType get_surface_type() const;

        // Optical CspElements setters.
        void set_aperture(const std::shared_ptr<Aperture> &aperture);
        void set_surface(const std::shared_ptr<Surface> &surface);
        void set_optics_front(const bool use_refraction, const float reflectivity,
                              const float transmissivity, const float slope_error, const float specularity_error,
                              const OpticalDistribution od);
        void set_optics_back(const bool use_refraction, const float reflectivity,
                             const float transmissivity, const float slope_error, const float specularity_error,
                             const OpticalDistribution od);
        void set_optics(
            const std::shared_ptr<SolTrace::Data::OpticalPropertySet>& optics);

        // return L2G rotation matrix
        Matrix33d get_rotation_matrix() const;
        void set_rotation_matrix(const Matrix33d& rotation_matrix);

        // return upper bounding box
        Vec3d get_upper_bounding_box() const;
        void set_upper_bounding_box(const Vec3d &upper);

        // return lower bounding box
        Vec3d get_lower_bounding_box() const;
        void set_lower_bounding_box(const Vec3d &lower);

        void set_bounding_box_local(const Vec3d &lower_local,
                                    const Vec3d &upper_local);

        // convert to device data available to GPU
        GeometryDataST toDeviceGeometryData() const override;
        MaterialData toDeviceMaterialDataFront() const override;
        MaterialData toDeviceMaterialDataBack() const override;

        // // we also need to implement the bounding box computation
        // // for a case like a rectangle aperture,
        // // once we have the origin, euler angles, rotatioin matrix
        // // and the aperture size, we can compute the bounding box
        // // this can be called when adding an element to the system
        // void compute_bounding_box();

        // Set id
        void set_id(const int32_t id) override;

    private:
        void set_optics(const bool is_front, const bool use_refraction, const float reflectivity,
                        const float transmissivity, const float slope_error, const float specularity_error,
                        const OpticalDistribution od);

        Vec3d m_origin;
        Vec3d m_aim_point;
        Matrix33d m_rotation_matrix;

        Vec3d m_upper_box_bound; // Global coordinates
        Vec3d m_lower_box_bound; // Global coordinates

        std::shared_ptr<Surface> m_surface;
        std::shared_ptr<Aperture> m_aperture;

        MaterialData m_optics_front;
        MaterialData m_optics_back;

        // Element id (from soltrace)
        int32_t m_id;
    };
}
