#include <cstdint>
#include <limits>
#include <vector>

#include "vec3d.h"
#include "soltrace_type.h"
#include "Surface.h"
#include "Aperture.h"
#include "utils/math_util.h"
#include "shaders/GeometryDataST.h"
#include "shaders/MaterialDataST.h"
#include "CspElement.h"
#include "soltrace_constants.h"

using namespace OptixCSP;

CspElementBase::CspElementBase()
{
}

CspElement::CspElement()
{
    m_origin = Vec3d(0.0, 0.0, 0.0);
    m_aim_point = Vec3d(0.0, 0.0, 1.0);
    m_rotation_matrix = Matrix33d();
    m_surface = nullptr;
    m_aperture = nullptr;
    m_id = kElementIdUnassigned;

    set_optics_front(false, 1.f, 0.f, 0.f, 0.f, OpticalDistribution::OPT_NONE);
    set_optics_back(false, 1.f, 0.f, 0.f, 0.f, OpticalDistribution::OPT_NONE);
}

// set and get origin
const Vec3d &CspElement::get_origin() const
{
    return m_origin;
}

void CspElement::set_origin(const Vec3d &o)
{
    m_origin = o;
}

const Vec3d &CspElement::get_aim_point() const
{
    return m_aim_point;
}

void CspElement::set_aim_point(const Vec3d &ap)
{
    m_aim_point = ap;
}

std::shared_ptr<Aperture> CspElement::get_aperture() const
{
    return m_aperture;
}

std::shared_ptr<Surface> CspElement::get_surface() const
{
    return m_surface;
}

ApertureType CspElement::get_aperture_type() const
{
    return m_aperture->get_aperture_type();
}

SurfaceType CspElement::get_surface_type() const
{
    return m_surface->get_surface_type();
}

// Optical CspElements setters.
void CspElement::set_aperture(const std::shared_ptr<Aperture> &aperture)
{
    m_aperture = aperture;
}
void CspElement::set_surface(const std::shared_ptr<Surface> &surface)
{
    m_surface = surface;
}

void CspElement::set_optics_front(const bool use_refraction, const float reflectivity,
                                  const float transmissivity, const float slope_error, const float specularity_error,
                                  const OpticalDistribution od)
{
    this->set_optics(true, use_refraction, reflectivity, transmissivity, slope_error, specularity_error, od);
}
void CspElement::set_optics_back(const bool use_refraction, const float reflectivity,
                                 const float transmissivity, const float slope_error, const float specularity_error,
                                 const OpticalDistribution od)
{
    this->set_optics(false, use_refraction, reflectivity, transmissivity, slope_error, specularity_error, od);
}

// return L2G rotation matrix
Matrix33d CspElement::get_rotation_matrix() const
{
    return m_rotation_matrix;
}

void CspElement::set_rotation_matrix(const Matrix33d &rotation_matrix)
{
    m_rotation_matrix = rotation_matrix;
}

// return upper bounding box
Vec3d CspElement::get_upper_bounding_box() const
{
    return m_upper_box_bound;
}

void CspElement::set_upper_bounding_box(const Vec3d &upper)
{
    m_upper_box_bound = upper;
    return;
}

// return lower bounding box
Vec3d CspElement::get_lower_bounding_box() const
{
    return m_lower_box_bound;
}

void CspElement::set_lower_bounding_box(const Vec3d &lower)
{
    m_lower_box_bound = lower;
    return;
}

void CspElement::set_bounding_box_local(const Vec3d &lower_local,
                                        const Vec3d &upper_local)
{
    // TODO: Functionality should be in simulation data. Perhaps with
    // a tighter bounding box. This version simply finds the global
    // bounding box for the local bounding box.

    Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix

    Vec3d c0 = lower_local;
    Vec3d c7 = upper_local;
    Vec3d c1(c0[0], c0[1], c7[2]);
    Vec3d c2(c0[0], c7[1], c0[2]);
    Vec3d c3(c7[0], c0[1], c0[2]);
    Vec3d c4(c0[0], c7[1], c7[2]);
    Vec3d c5(c7[0], c0[1], c7[2]);
    Vec3d c6(c7[0], c7[1], c0[2]);

    Vec3d g0 = rotation_matrix * c0 + m_origin;
    Vec3d g1 = rotation_matrix * c1 + m_origin;
    Vec3d g2 = rotation_matrix * c2 + m_origin;
    Vec3d g3 = rotation_matrix * c3 + m_origin;
    Vec3d g4 = rotation_matrix * c4 + m_origin;
    Vec3d g5 = rotation_matrix * c5 + m_origin;
    Vec3d g6 = rotation_matrix * c6 + m_origin;
    Vec3d g7 = rotation_matrix * c7 + m_origin;

    // go through the corners and find the min and max x, y, z
    std::vector<Vec3d> corners = {g0, g1, g2, g3,
                                  g4, g5, g6, g7};

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();

    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();

    for (auto &corner : corners)
    {
        min_x = fmin(min_x, corner[0]);
        min_y = fmin(min_y, corner[1]);
        min_z = fmin(min_z, corner[2]);

        max_x = fmax(max_x, corner[0]);
        max_y = fmax(max_y, corner[1]);
        max_z = fmax(max_z, corner[2]);
    }

    // set the lower and upper bounds
    m_lower_box_bound[0] = min_x;
    m_lower_box_bound[1] = min_y;
    m_lower_box_bound[2] = min_z;

    m_upper_box_bound[0] = max_x;
    m_upper_box_bound[1] = max_y;
    m_upper_box_bound[2] = max_z;

    return;
}

GeometryDataST CspElement::toDeviceGeometryData() const
{

    GeometryDataST geometry_data;

    SurfaceType surface_type = m_surface->get_surface_type();
    ApertureType aperture_type = m_aperture->get_aperture_type();

    if (aperture_type == ApertureType::RECTANGLE)
    {

        double width = m_aperture->get_width();
        double height = m_aperture->get_height();
        auto rect_aperture = std::dynamic_pointer_cast<ApertureRectangle>(m_aperture);
        double x_coord = rect_aperture->get_x_coord();
        double y_coord = rect_aperture->get_y_coord();

        Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix

        Vec3d v1 = rotation_matrix.get_x_basis();
        Vec3d v2 = rotation_matrix.get_y_basis();

        if (surface_type == SurfaceType::FLAT)
        {
            Vec3d local_center(x_coord + 0.5 * width, y_coord + 0.5 * height, 0);
            Vec3d global_center = rotation_matrix * local_center + m_origin;
            GeometryDataST::Rectangle_Flat heliostat(OptixCSP::toFloat3(global_center), OptixCSP::toFloat3(v1), OptixCSP::toFloat3(v2), (float)width, (float)height);
            geometry_data.setRectangle_Flat(heliostat);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            Vec3d edge_x = v1 * (float)(width);
            Vec3d edge_y = v2 * (float)(height);

            // Lower left corner
            Vec3d local_anchor(x_coord, y_coord, 0.0);
            // float3 anchor = OptixCSP::toFloat3(m_origin - v1 * 0.5 - v2 * 0.5);
            Vec3d global_anchor = rotation_matrix * local_anchor + m_origin;

            GeometryDataST::Rectangle_Parabolic heliostat(OptixCSP::toFloat3(edge_x),
                                                          OptixCSP::toFloat3(edge_y),
                                                          OptixCSP::toFloat3(global_anchor),
                                                          (float)m_surface->get_curvature_1(),
                                                          (float)m_surface->get_curvature_2());
            geometry_data.setRectangleParabolic(heliostat);
        }

        if (surface_type == SurfaceType::CYLINDER)
        {
            float radius = static_cast<float>(width) / 2.0f;
            float half_height = static_cast<float>(height) / 2.0f;

            float3 center = OptixCSP::toFloat3(m_origin);
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            float3 base_x = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 base_z = OptixCSP::toFloat3(rotation_matrix.get_z_basis());

            GeometryDataST::Cylinder_Y heliostat(center, radius, half_height, base_x, base_z);
            geometry_data.setCylinder_Y(heliostat);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            float3 vx = OptixCSP::toFloat3(v1);
            float3 vy = OptixCSP::toFloat3(v2);
            float3 o  = OptixCSP::toFloat3(m_origin);
            float R = 1.0f / (float)m_surface->get_curvature_1();
            GeometryDataST::Rectangle_Spherical rs(o, vx, vy, R,
                                                   (float)width, (float)height,
                                                   (float)x_coord, (float)y_coord);
            geometry_data.setRectangle_Spherical(rs);
        }
    }

    if (aperture_type == ApertureType::TRIANGLE)
    {
        Vec3d v1, v2, v3;
        // first cast to ApertureTriangle type
        ApertureTriangle tri = static_cast<ApertureTriangle &>(*m_aperture);

        v1 = tri.get_v0();
        v2 = tri.get_v1();
        v3 = tri.get_v2();

        if (surface_type == SurfaceType::FLAT)
        {
            // given the origin and rotation, compute global coordinates of the triangle vertices
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            Vec3d v1_global = rotation_matrix * v1 + m_origin;
            Vec3d v2_global = rotation_matrix * v2 + m_origin;
            Vec3d v3_global = rotation_matrix * v3 + m_origin;

            GeometryDataST::Triangle_Flat heliostat(OptixCSP::toFloat3(v1_global),
                                                    OptixCSP::toFloat3(v2_global),
                                                    OptixCSP::toFloat3(v3_global));
            geometry_data.setTriangle_Flat(heliostat);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());

            float cx = (float)(m_surface->get_curvature_1());
            float cy = (float)(m_surface->get_curvature_2());
            float3 o = OptixCSP::toFloat3(m_origin);

            // Triangle vertices are in local element XY frame (z=0).
            // The 2D local x,y coords for the aperture test equal the
            // Vec3d x and y components directly (since rotation is orthonormal).
            float2 lv0 = make_float2((float)v1[0], (float)v1[1]);
            float2 lv1 = make_float2((float)v2[0], (float)v2[1]);
            float2 lv2 = make_float2((float)v3[0], (float)v3[1]);

            GeometryDataST::Triangle_Parabolic trip(o, vx, vy, cx, cy, lv0, lv1, lv2);
            geometry_data.setTriangle_Parabolic(trip);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            Matrix33d rotation_matrix = get_rotation_matrix();
            float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());
            float R   = 1.0f / (float)(m_surface->get_curvature_1());
            float3 o  = OptixCSP::toFloat3(m_origin);

            float2 lv0 = make_float2((float)v1[0], (float)v1[1]);
            float2 lv1 = make_float2((float)v2[0], (float)v2[1]);
            float2 lv2 = make_float2((float)v3[0], (float)v3[1]);

            GeometryDataST::Triangle_Spherical tris(o, vx, vy, R, lv0, lv1, lv2);
            geometry_data.setTriangle_Spherical(tris);
        }
    }

    if (aperture_type == ApertureType::QUADRILATERAL)
    {
        ApertureQuadrilateral quad = static_cast<ApertureQuadrilateral &>(*m_aperture);

        Vec3d p1, p2, p3, p4;

        p1 = quad.get_p0();
        p2 = quad.get_p1();
        p3 = quad.get_p2();
        p4 = quad.get_p3();

        if (surface_type == SurfaceType::FLAT)
        {
            // given the origin and rotation, compute global coordinates of the triangle vertices
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            Vec3d p1_global = rotation_matrix * p1 + m_origin;
            Vec3d p2_global = rotation_matrix * p2 + m_origin;
            Vec3d p3_global = rotation_matrix * p3 + m_origin;
            Vec3d p4_global = rotation_matrix * p4 + m_origin;

            GeometryDataST::Quadrilateral_Flat heliostat(OptixCSP::toFloat3(p1_global),
                                                         OptixCSP::toFloat3(p2_global),
                                                         OptixCSP::toFloat3(p3_global),
                                                         OptixCSP::toFloat3(p4_global));
            geometry_data.setQuadrilateral_Flat(heliostat);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());

            float cx = (float)(m_surface->get_curvature_1());
            float cy = (float)(m_surface->get_curvature_2());
            float3 o = OptixCSP::toFloat3(m_origin);

            GeometryDataST::Quadrilateral_Parabolic heliostat(
                o, vx, vy, cx, cy,
                make_float2(p1[0], p1[1]), make_float2(p2[0], p2[1]), 
                make_float2(p3[0], p3[1]), make_float2(p4[0], p4[1]));
            geometry_data.setQuadrilateral_Parabolic(heliostat);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            Matrix33d rotation_matrix = get_rotation_matrix();
            float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());
            float R   = 1.0f / (float)(m_surface->get_curvature_1());
            float3 o  = OptixCSP::toFloat3(m_origin);

            GeometryDataST::Quadrilateral_Spherical qus(
                o, vx, vy, R,
                make_float2(p1[0], p1[1]), make_float2(p2[0], p2[1]),
                make_float2(p3[0], p3[1]), make_float2(p4[0], p4[1]));
            geometry_data.setQuadrilateral_Spherical(qus);
        }
    }

    if (aperture_type == ApertureType::CIRCLE)
    {
        ApertureCircle circ = static_cast<ApertureCircle &>(*m_aperture);
        float r = circ.get_radius();
        float3 o = OptixCSP::toFloat3(m_origin);
        float3 n = normalize(OptixCSP::toFloat3(m_aim_point - m_origin));

        if (surface_type == SurfaceType::FLAT)
        {
            GeometryDataST::Circle_Flat heliostat(o, n, r);
            geometry_data.setCircle_Flat(heliostat);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
            float3 v1 = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 v2 = OptixCSP::toFloat3(rotation_matrix.get_y_basis());
            float cx = (float)(m_surface->get_curvature_1());
            float cy = (float)(m_surface->get_curvature_2());
            GeometryDataST::Circle_Parabolic heliostat(o, v1, v2, cx, cy, r);
            geometry_data.setCircle_Parabolic(heliostat);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            Matrix33d rotation_matrix = get_rotation_matrix();
            float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
            float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());
            float R   = 1.0f / (float)(m_surface->get_curvature_1());
            GeometryDataST::Circle_Spherical cs(o, vx, vy, R, r);
            geometry_data.setCircle_Spherical(cs);
        }
    }

    if (aperture_type == ApertureType::HEXAGON)
    {
        Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
        float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
        float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());

        ApertureHexagon hex = static_cast<ApertureHexagon &>(*m_aperture);
        float s = hex.get_side_length();
        float3 o = OptixCSP::toFloat3(m_origin);
        float3 n = normalize(OptixCSP::toFloat3(m_aim_point - m_origin));

        if (surface_type == SurfaceType::FLAT)
        {
            GeometryDataST::Hexagon_Flat hex(o, n, vx, vy, s);
            geometry_data.setHexagon_Flat(hex);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            float cx = (float)(m_surface->get_curvature_1());
            float cy = (float)(m_surface->get_curvature_2());
            GeometryDataST::Hexagon_Parabolic hex(o, vx, vy, cx, cy, s);
            geometry_data.setHexagon_Parabolic(hex);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            float R = 1.0f / (float)(m_surface->get_curvature_1());
            GeometryDataST::Hexagon_Spherical hs(o, vx, vy, R, s);
            geometry_data.setHexagon_Spherical(hs);
        }
    }

    if (aperture_type == ApertureType::ANNULUS)
    {
        Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix
        float3 vx = OptixCSP::toFloat3(rotation_matrix.get_x_basis());
        float3 vy = OptixCSP::toFloat3(rotation_matrix.get_y_basis());

        ApertureAnnulus anf = static_cast<ApertureAnnulus &>(*m_aperture);
        float radius_in = anf.get_radius_inner();
        float radius_out = anf.get_radius_outer();
        float arc = anf.get_arc();
        float3 o = OptixCSP::toFloat3(m_origin);
        float3 n = normalize(OptixCSP::toFloat3(m_aim_point - m_origin));

        if (surface_type == SurfaceType::FLAT)
        {
            GeometryDataST::Annulus_Flat anf(o, n, vx, vy,
                                             radius_in, radius_out, arc);
            geometry_data.setAnnulus_Flat(anf);
        }

        if (surface_type == SurfaceType::PARABOLIC)
        {
            float cx = (float)(m_surface->get_curvature_1());
            float cy = (float)(m_surface->get_curvature_2());
            GeometryDataST::Annulus_Parabolic anp(o, vx, vy, cx, cy,
                                                  radius_in, radius_out, arc);
            geometry_data.setAnnulus_Parabolic(anp);
        }

        if (surface_type == SurfaceType::SPHERICAL)
        {
            float R = 1.0f / (float)(m_surface->get_curvature_1());
            GeometryDataST::Annulus_Spherical as(o, vx, vy, R,
                                                 radius_in, radius_out, arc);
            geometry_data.setAnnulus_Spherical(as);
        }
    }

    geometry_data.id = this->m_id;

    return geometry_data;
}

MaterialData CspElement::toDeviceMaterialDataFront() const
{
    return this->m_optics_front;
}

MaterialData CspElement::toDeviceMaterialDataBack() const
{
    return this->m_optics_back;
}

// // we also need to implement the bounding box computation
// // for a case like a rectangle aperture,
// // once we have the origin, euler angles, rotation matrix
// // and the aperture size, we can compute the bounding box
// // this can be called when adding an element to the system
// void CspElement::compute_bounding_box()
// {
//     // this can also be called while "initializing" the element
//     // get the rotation matrix first
//     Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix

//     // now check the type of the aperture
//     ApertureType aperture_type = m_aperture->get_aperture_type();
//     SurfaceType surface_type = m_surface->get_surface_type();

//     if (aperture_type == ApertureType::RECTANGLE && surface_type != SurfaceType::CYLINDER)
//     {
//         // get the width and height of the aperture
//         double width = m_aperture->get_width();
//         double height = m_aperture->get_height();

//         // compute the four corners of the rectangle locally
//         Vec3d corner1 = Vec3d(-width / 2, -height / 2, 0.0);
//         Vec3d corner2 = Vec3d(width / 2, -height / 2, 0.0);
//         Vec3d corner3 = Vec3d(width / 2, height / 2, 0.0);
//         Vec3d corner4 = Vec3d(-width / 2, height / 2, 0.0);

//         // transform the corners to the global frame
//         Vec3d corner1_global = rotation_matrix * corner1 + m_origin;
//         Vec3d corner2_global = rotation_matrix * corner2 + m_origin;
//         Vec3d corner3_global = rotation_matrix * corner3 + m_origin;
//         Vec3d corner4_global = rotation_matrix * corner4 + m_origin;

//         // Extended z axis box slightly
//         double epsilon = 1e-1;

//         // now update the bounding box, need to find the min and max x, y, z
//         m_lower_box_bound[0] = fmin(fmin(corner1_global[0], corner2_global[0]), fmin(corner3_global[0], corner4_global[0]));
//         m_lower_box_bound[1] = fmin(fmin(corner1_global[1], corner2_global[1]), fmin(corner3_global[1], corner4_global[1]));
//         m_lower_box_bound[2] = fmin(fmin(corner1_global[2], corner2_global[2]), fmin(corner3_global[2], corner4_global[2])) - epsilon;

//         m_upper_box_bound[0] = fmax(fmax(corner1_global[0], corner2_global[0]), fmax(corner3_global[0], corner4_global[0]));
//         m_upper_box_bound[1] = fmax(fmax(corner1_global[1], corner2_global[1]), fmax(corner3_global[1], corner4_global[1]));
//         m_upper_box_bound[2] = fmax(fmax(corner1_global[2], corner2_global[2]), fmax(corner3_global[2], corner4_global[2])) + epsilon;
//     }

//     // slightly different for the cylinder, we want to know the radius and half height
//     if (surface_type == SurfaceType::CYLINDER)
//     {
//         // get the radius and full height of the cylinder
//         double width = m_aperture->get_width();
//         double height = m_aperture->get_height();

//         // compute 8 corners of the cyliinder box locally
//         Vec3d corner1 = Vec3d(-width / 2, -height / 2, -width / 2);
//         Vec3d corner2 = Vec3d(width / 2, -height / 2, -width / 2);
//         Vec3d corner3 = Vec3d(-width / 2, height / 2, -width / 2);
//         Vec3d corner4 = Vec3d(width / 2, height / 2, -width / 2);
//         Vec3d corner5 = Vec3d(-width / 2, -height / 2, width / 2);
//         Vec3d corner6 = Vec3d(width / 2, -height / 2, width / 2);
//         Vec3d corner7 = Vec3d(-width / 2, height / 2, width / 2);
//         Vec3d corner8 = Vec3d(width / 2, height / 2, width / 2);

//         // get the rotation matrix
//         Matrix33d rotation_matrix = get_rotation_matrix(); // L2G rotation matrix

//         // transform the corners to the global frame
//         Vec3d corner1_global = rotation_matrix * corner1 + m_origin;
//         Vec3d corner2_global = rotation_matrix * corner2 + m_origin;
//         Vec3d corner3_global = rotation_matrix * corner3 + m_origin;
//         Vec3d corner4_global = rotation_matrix * corner4 + m_origin;
//         Vec3d corner5_global = rotation_matrix * corner5 + m_origin;
//         Vec3d corner6_global = rotation_matrix * corner6 + m_origin;
//         Vec3d corner7_global = rotation_matrix * corner7 + m_origin;
//         Vec3d corner8_global = rotation_matrix * corner8 + m_origin;

//         // go through the corners and find the min and max x, y, z
//         std::vector<Vec3d> corners = {corner1_global, corner2_global, corner3_global, corner4_global,
//                                       corner5_global, corner6_global, corner7_global, corner8_global};

//         double min_x = std::numeric_limits<double>::max();
//         double min_y = std::numeric_limits<double>::max();
//         double min_z = std::numeric_limits<double>::max();

//         double max_x = std::numeric_limits<double>::lowest();
//         double max_y = std::numeric_limits<double>::lowest();
//         double max_z = std::numeric_limits<double>::lowest();

//         for (auto &corner : corners)
//         {
//             min_x = fmin(min_x, corner[0]);
//             min_y = fmin(min_y, corner[1]);
//             min_z = fmin(min_z, corner[2]);

//             max_x = fmax(max_x, corner[0]);
//             max_y = fmax(max_y, corner[1]);
//             max_z = fmax(max_z, corner[2]);
//         }

//         // set the lower and upper bounds
//         m_lower_box_bound[0] = min_x;
//         m_lower_box_bound[1] = min_y;
//         m_lower_box_bound[2] = min_z;

//         m_upper_box_bound[0] = max_x;
//         m_upper_box_bound[1] = max_y;
//         m_upper_box_bound[2] = max_z;
//     }

//     // bounding box for triangle aperture
//     if (aperture_type == ApertureType::TRIANGLE)
//     {
//         // get the three vertices of the triangle aperture in local coordinates
//         // first cast to ApertureTriangle type
//         ApertureTriangle tri = static_cast<ApertureTriangle &>(*m_aperture);

//         Vec3d v1 = tri.get_v0();
//         Vec3d v2 = tri.get_v1();
//         Vec3d v3 = tri.get_v2();
//         // transform the vertices to the global frame
//         Vec3d v1_global = rotation_matrix * v1 + m_origin;
//         Vec3d v2_global = rotation_matrix * v2 + m_origin;
//         Vec3d v3_global = rotation_matrix * v3 + m_origin;
//         // now update the bounding box, need to find the min and max x, y, z
//         m_lower_box_bound[0] = fmin(fmin(v1_global[0], v2_global[0]), v3_global[0]);
//         m_lower_box_bound[1] = fmin(fmin(v1_global[1], v2_global[1]), v3_global[1]);
//         m_lower_box_bound[2] = fmin(fmin(v1_global[2], v2_global[2]), v3_global[2]);
//         m_upper_box_bound[0] = fmax(fmax(v1_global[0], v2_global[0]), v3_global[0]);
//         m_upper_box_bound[1] = fmax(fmax(v1_global[1], v2_global[1]), v3_global[1]);
//         m_upper_box_bound[2] = fmax(fmax(v1_global[2], v2_global[2]), v3_global[2]);
//     }
// }

void CspElement::set_id(const int32_t id)
{
    this->m_id = id;
}

void CspElement::set_optics(const bool is_front, const bool use_refraction, const float reflectivity,
                            const float transmissivity, const float slope_error, const float specularity_error,
                            const OpticalDistribution od)
{
    auto &md = is_front ? this->m_optics_front : this->m_optics_back;
    md.use_refraction = use_refraction;
    md.reflectivity = reflectivity;
    md.transmissivity = transmissivity;
    md.slope_error = slope_error;
    md.specularity_error = specularity_error;
    md.optical_dist = od;
}
