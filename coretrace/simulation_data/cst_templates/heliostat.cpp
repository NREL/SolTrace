
#include "heliostat.hpp"

#include <sstream>
#include <stdexcept>

#include "utilities.hpp"

namespace SolTrace::Data
{

    Heliostat::Heliostat()
        : CompositeElement(),
          initialized(false),
          aperture_size_x(-1.0),
          aperture_size_y(-1.0),
          focal_length_x(-1.0),
          focal_length_y(-1.0),
          gap_x(-1.0),
          gap_y(-1.0),
          num_panels_x(-1),
          num_panels_y(-1),
          canting_method(UNSET),
          onaxis_canting_distance(-1.0),
          offaxis_canting_sun_position_azimuth(-1.0),
          offaxis_canting_sun_position_zenith(-1.0),
          heliostat_area(-1.0),
          tracking_azimuth(-1.0),
          tracking_elevation(-1.0),
          target_set(false)
    {
        this->elevation_axis.set_values(1.0, 0.0, 0.0);
        this->sun_position.set_values(0.0, 0.0, 1.0);
        this->optics_mirror.set_ideal_reflection();
        return;
    }

    Heliostat::~Heliostat()
    {
        this->facets.clear();
        return;
    }

    void Heliostat::create_geometry()
    {
        if (this->initialized)
            return;

        this->enforce_user_fields_set();

        this->facets.clear();
        this->clear();

        double panel_len_x = this->aperture_size_x -
                             this->gap_x * (this->num_panels_x - 1);
        panel_len_x /= this->num_panels_x;
        double panel_len_y = this->aperture_size_y -
                             this->gap_y * (this->num_panels_y - 1);
        panel_len_y /= this->num_panels_y;

        if (this->canting_method == OFF_AXIS)
        {
            // TODO: Do stuff here...
        }

        this->heliostat_area = 0.0;
        double panel_x = -0.5 * (this->aperture_size_x - panel_len_x);
        double panel_y = -0.5 * (this->aperture_size_y - panel_len_y);

        single_element_ptr elem;
        Vector3d origin;
        Vector3d aim;
        aperture_ptr ap;
        surface_ptr surf;
        element_id sts;
        double c, z;

        for (int j = 0; j < this->num_panels_y; ++j)
        {
            for (int i = 0; i < this->num_panels_x; ++i)
            {
                elem = make_element<SingleElement>();

                if (this->canting_method == NONE)
                {
                    origin.set_values(panel_x, panel_y, 0.0);
                    aim.set_values(panel_x, panel_y, 1000.0);
                }
                else if (this->canting_method == ON_AXIS)
                {
                    c = 0.5 / (this->onaxis_canting_distance);
                    z = 0.5 * c * (panel_x * panel_x + panel_y * panel_y);
                    origin.set_values(panel_x, panel_y, z);
                    aim.set_values(0.0, 0.0, 2.0 * this->onaxis_canting_distance);
                }
                else if (this->canting_method == OFF_AXIS)
                {
                    origin.set_values(panel_x, panel_y, 0.0);
                    // TODO: Set aim vector values
                    aim.set_values(0.0, 0.0, 1.0);
                    throw std::runtime_error("OFF_AXIS is not yet implemented");
                }
                else if (this->canting_method == UNSET)
                {
                    throw std::runtime_error("Canting method unset");
                }
                else
                {
                    throw std::runtime_error("Unknown canting method");
                }

                ap = make_aperture<Rectangle>(panel_len_x, panel_len_y);
                elem->set_aperture(ap);

                if (this->focal_length_x <= 0.0 && this->focal_length_y <= 0.0)
                {
                    surf = make_surface<Flat>();
                }
                else
                {
                    surf = make_surface<Parabola>(this->focal_length_x,
                                                  this->focal_length_y);
                }
                elem->set_surface(surf);

                elem->set_reference_frame_geometry(origin,
                                                   aim,
                                                   0.0);

                // TODO: Make back optical properties accessible to user
                elem->set_front_optical_properties(this->optics_mirror);
                elem->set_back_optical_properties(this->optics_mirror);
                elem->enable();

                this->heliostat_area += panel_len_x * panel_len_y;
                this->facets.push_back(elem);
                sts = this->add_element(elem);
                if (!Element::is_success(sts))
                {
                    throw std::runtime_error("Failed to add element");
                }

                panel_x += panel_len_x + this->gap_x;
            }
            panel_x = -0.5 * (this->aperture_size_x - panel_len_x);
            panel_y += panel_len_y + this->gap_y;
        }

        this->initialized = true;

        return;
    }

    void Heliostat::update_geometry(double azimuth,
                                    double elevation)
    {

        if (elevation < 0.0 || elevation > 90.0)
        {
            std::stringstream ss;
            ss << "Heliostat::update_geometry: Invalid elevation ("
               << elevation << "). Elevation must lie between 0 and 90 degrees.";
            throw std::invalid_argument(ss.str());
        }

        if (azimuth < -180.0 || azimuth > 180.0)
        {
            std::stringstream ss;
            ss << "Heliostat::update_geometry: Invalid azimuth ("
               << elevation << "). Azimuth must lie between -180 and 180 degrees.";
            throw std::invalid_argument(ss.str());
        }

        if (!this->initialized)
        {
            throw std::runtime_error(
                "Heliostat::update_geometry: Must call create_geometry prior "
                "to update_geometry");
        }

        this->coordinates_initialized = false;

        this->tracking_azimuth = azimuth;
        this->tracking_elevation = elevation;

        sun_position_vector_degrees(this->sun_position, azimuth, elevation);
        Vector3d target_dir;
        vector_add(1.0, this->target_pos,
                   -1.0, this->get_origin_global(),
                   target_dir);

        Vector3d aim_vector;
        this->sun_position.make_unit();
        target_dir.make_unit();
        vector_add(1.0, target_dir, 1.0, this->sun_position, aim_vector);
        aim_vector.make_unit();
        this->convert_global_to_reference(this->aim, aim_vector);
        aim_vector = this->aim;
        vector_add(1.0, this->get_origin_ref(),
                   1000.0, this->aim);

        // Project into xy-plane
        aim_vector[2] = 0.0;
        double theta = acos(aim_vector[0] / vector_norm(aim_vector));
        this->set_zrot_radians(theta);

        // std::cout << "Origin: " << this->origin
        //           << "\nAim Point: " << this->aim
        //           << "\nZ Rot: " << this->zrot
        //           << "\nAim Vector: " << aim_vector
        //           << "\nTarget: " << this->target_pos
        //           << "\nTarget Vector: " << target_dir
        //           << std::endl;

        this->compute_coordinate_rotations();

        return;
    }

    void Heliostat::set_aperture_size(double size_x,
                                      double size_y)
    {
        if (size_x <= 0.0 || size_y <= 0.0)
        {
            std::stringstream ss;
            ss << "Heliostat::set_aperture_size: Invalid aperture dimensions. "
               << "size_x (" << size_x << ") and size_y (" << size_y
               << ") must be positive.";
            throw std::invalid_argument(ss.str());
        }

        this->initialized = false;
        this->aperture_size_x = size_x;
        this->aperture_size_y = size_y;

        return;
    }

    void Heliostat::set_focal_length(double flen)
    {
        if (flen < 0.0)
        {
            std::stringstream ss;
            ss << "Heliostat::set_focal_length: Invalid focal length ("
               << flen << "). Must be non-negative.";
            throw std::invalid_argument(ss.str());
        }

        this->initialized = false;
        this->focal_length_x = flen;
        this->focal_length_y = flen;

        return;
    }

    void Heliostat::set_focal_length(double fx, double fy)
    {
        if (fx < 0.0 || fy < 0.0)
        {
            std::stringstream ss;
            ss << "Heliostat::set_focal_length: Invalid focal lengths. "
               << "fx (" << fx << ") and fy (" << fy
               << ") must be non-negative.";
            throw std::invalid_argument(ss.str());
        }

        this->initialized = false;
        this->focal_length_x = fx;
        this->focal_length_y = fy;

        return;
    }

    void Heliostat::set_gaps(double gap_x,
                             double gap_y)
    {
        if (gap_x < 0.0 || gap_y < 0.0)
        {
            std::stringstream ss;
            ss << "Heliostat::set_gaps: Invalid gap dimensions. "
               << "gap_x (" << gap_x << ") and gap_y (" << gap_y
               << ") must be non-negative.";
            throw std::invalid_argument(ss.str());
        }

        this->initialized = false;
        this->gap_x = gap_x;
        this->gap_y = gap_y;

        return;
    }

    void Heliostat::set_number_panels(uint_fast64_t num_x,
                                      uint_fast64_t num_y)
    {
        if (num_x < 1 || num_y < 1)
        {
            std::stringstream ss;
            ss << "Heliostat::set_number_panels: Invalid panel count. "
               << "num_x (" << num_x << ") and num_y (" << num_y
               << ") must be at least 1.";
            throw std::invalid_argument(ss.str());
        }

        this->initialized = false;
        this->num_panels_x = num_x;
        this->num_panels_y = num_y;

        return;
    }

    void Heliostat::set_optics(const OpticalProperties &optics)
    {
        this->optics_mirror = optics;
        // TODO: Need to update the subelements!
        return;
    }

    void Heliostat::set_canting(CantingType ct, double val1, double val2)
    {
        this->initialized = false;
        this->canting_method = ct;

        if (ct == NONE)
        {
            ; // Intentional no-op
        }
        else if (ct == OFF_AXIS)
        {
            if (val1 < 0.0 || val1 > 360.0)
            {
                std::stringstream ss;
                ss << "Heliostat::set_canting: Invalid off-axis canting "
                   << "sun azimuth angle (" << val1
                   << "). Must be between 0 and 360 degrees.";
                throw std::invalid_argument(ss.str());
            }
            if (val2 < 0.0 || val2 > 90.0)
            {
                std::stringstream ss;
                ss << "Heliostat::set_canting: Invalid off-axis canting "
                   << "sun zenith angle (" << val2
                   << "). Must be between 0 and 90 degrees.";
                throw std::invalid_argument(ss.str());
            }
            this->offaxis_canting_sun_position_azimuth = val1;
            this->offaxis_canting_sun_position_zenith = val2;
            this->onaxis_canting_distance = -1.0;
        }
        else if (ct == ON_AXIS)
        {
            if (val1 <= 0.0)
            {
                std::stringstream ss;
                ss << "Heliostat::set_canting: Invalid on-axis canting "
                   << "distance (" << val1 << "). Must be positive.";
                throw std::invalid_argument(ss.str());
            }
            this->onaxis_canting_distance = val1;
            this->offaxis_canting_sun_position_azimuth = -1.0;
            this->offaxis_canting_sun_position_zenith = -1.0;
        }
        else
        {
            std::stringstream ss;
            ss << "Heliostat::set_canting: Unrecognized canting method ("
               << static_cast<int>(ct) << ").";
            throw std::invalid_argument(ss.str());
        }

        return;
    }

    void Heliostat::set_target_position(const Vector3d &pos)
    {
        this->target_pos = pos;
        this->target_set = true;
        return;
    }

    void Heliostat::set_tracking_limits(double az_lower, double az_upper,
                                        double el_lower, double el_upper)
    {
        // TODO: Implement this
        return;
    }

    void Heliostat::enforce_user_fields_set() const
    {
        if (this->initialized)
        {
            CompositeElement::enforce_user_fields_set();
        }

        // Validate that all required parameters have been set
        if (this->aperture_size_x <= 0.0 || this->aperture_size_y <= 0.0)
        {
            throw std::invalid_argument(
                "Heliostat: Aperture size must be set "
                "before creating geometry.");
        }

        if (this->num_panels_x <= 0 || this->num_panels_y <= 0)
        {
            throw std::invalid_argument(
                "Heliostat: Number of panels must be set "
                "before creating geometry.");
        }

        if (this->gap_x < 0.0 || this->gap_y < 0.0)
        {
            throw std::invalid_argument(
                "Heliostat: Panel gap must be set "
                "before creating geometry.");
        }

        if (this->canting_method == UNSET)
        {
            throw std::invalid_argument(
                "Heliostat: Canting method must be set "
                "before creating geometry.");
        }

        if (!this->target_set)
        {
            throw std::invalid_argument(
                "Heliostat: Target position must be set "
                "before creating geometry.");
        }

        return;
    }

} // namespace SolTrace::Data
