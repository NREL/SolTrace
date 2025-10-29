
#include "aperture.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

// #include <iostream>

#include "constants.hpp"

namespace SolTrace::Data {

aperture_ptr Aperture::make_aperture_from_type(ApertureType type,
                                               const std::vector<double> &args)
{
    // std::cout << "Type: " << type
    //           << " NArgs: " << args.size()
    //           << std::endl;
    switch (type)
    {
    case ApertureType::ANNULUS:
        if (args.size() < 3)
            break;
        return make_aperture<Annulus>(args[0], args[1], args[2]);
    case ApertureType::CIRCLE:
        if (args.size() < 1)
            break;
        return make_aperture<Circle>(args[0]);
    case ApertureType::HEXAGON:
        if (args.size() < 1)
            break;
        return make_aperture<Hexagon>(args[0]);
    case ApertureType::RECTANGLE:
        if (args.size() < 2)
            break;
        return make_aperture<Rectangle>(args[0], args[1]); // This is assuming centered around the origin
    case ApertureType::EQUILATERAL_TRIANGLE:
        if (args.size() < 1)
            break;
        return make_aperture<EqualateralTriangle>(args[0]);
    case ApertureType::SINGLE_AXIS_CURVATURE_SECTION:
        if (args.size() < 3)
            break;
        return make_aperture<Rectangle>(
            args[1] - args[0], args[2], 
            -0.5 * (args[1] - args[0]), -0.5 * args[2]);
        // return make_aperture<Rectangle>(args[1] - args[0], args[2]);
    case ApertureType::IRREGULAR_TRIANGLE:
        if (args.size() < 6)
            break;
        return make_aperture<IrregularTriangle>(
            args[0], args[1], args[2], args[3], args[4], args[5]);
    case ApertureType::IRREGULAR_QUADRILATERAL:
        if (args.size() < 8)
            break;
        return make_aperture<IrregularQuadrilateral>(
            args[0], args[1], args[2], args[3],
            args[4], args[5], args[6], args[7]);
    default:
        // TODO handle error
        // Unsupported case
        return nullptr;
        // break;
    }

    // TODO handle error
    // Wrong number of arguments

    return nullptr;
    // return aperture_ptr();
}

double Annulus::aperture_area() const
{
    // TODO: input.cpp on line 219 uses the formula
    //    elm->ParameterC*(ACOSM1O180)*(elm->ParameterB - elm->ParameterA);
    //    = \theta * (r - R)
    // This seems to be wrong...
    double R = this->outer_radius;
    double r = this->inner_radius;
    // Convert to radians
    double arc = this->arc_angle * D2R;
    return 0.5 * arc * (R * R - r * r);
}

double Annulus::diameter_circumscribed_circle() const
{
    return 2.0 * this->outer_radius;
}

bool Annulus::is_in(double x, double y) const
{
    double r = sqrt(x * x + y * y);
    bool inside = false;
    if (this->inner_radius <= r &&
        r <= this->outer_radius)
    {
        double theta = atan2(y, x);
        // Arc is split across x-axis, hence the 0.5
        double arc = 0.5 * this->arc_angle * D2R;
        inside = (-arc <= theta && theta <= arc);
    }
    return inside;
}

aperture_ptr Annulus::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<Annulus>(*this);
}

double Circle::aperture_area() const
{
    return 0.25 * PI * this->diameter * this->diameter;
}

double Circle::diameter_circumscribed_circle() const
{
    return this->diameter;
}

bool Circle::is_in(double x, double y) const
{
    double r = sqrt(x * x + y * y);
    return r <= this->radius_circumscribed_circle();
}

aperture_ptr Circle::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<Circle>(*this);
}

double EqualateralTriangle::aperture_area() const
{
    double r = 0.5 * this->circumscribe_diameter;
    return 0.75 * sqrt(3) * r * r;
}

double EqualateralTriangle::diameter_circumscribed_circle() const
{
    return this->circumscribe_diameter;
}

bool EqualateralTriangle::is_in(double x, double y) const
{
    double r = sqrt(x * x + y * y);
    double ro = this->radius_circumscribed_circle();
    if (r > ro)
        return false;

    double ri = 0.5 * ro;
    if (r <= ri)
        return true;

    double y0;
    // double a = ro / sqrt(3.0) = 2 * ri / sqrt(3.0);
    if (0.0 <= x && x <= ro)
    {
        // y0 = -sqrt(3.0) * (x - a);
        y0 = ro - sqrt(3.0) * x;
        return (-ri <= y && y <= y0);
    }
    else if (-ro <= x && x < 0.0)
    {
        // y0 = sqrt(3.0) * (x + a);
        y0 = sqrt(3.0) * x + ro;
        return (-ri <= y && y <= y0);
    }

    return false;
}

aperture_ptr EqualateralTriangle::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<EqualateralTriangle>(*this);
}

double Hexagon::aperture_area() const
{
    // TODO: input.cpp on line 210 uses the formula
    //    5*sqr(elm->ParameterA/2.0)*cos(30.0*(ACOSM1O180))*sin(30.0*(ACOSM1O180));
    //    = 5*(d/2)^2*cos(pi/6)*sin(pi/6)
    //    = 5*(d/2)^2*sqrt(3)/2*1/2
    //    = 5*sqrt(3)/4 * (d/2)^2
    //    = 1.25*sqrt(3) * (d/2)^2
    // This seems to be wrong...
    double r = 0.5 * this->circumscribe_diameter;
    return 1.5 * sqrt(3) * r * r;
}

double Hexagon::diameter_circumscribed_circle() const
{
    return circumscribe_diameter;
}

bool Hexagon::is_in(double x, double y) const
{
    double r = sqrt(x * x + y * y);
    double ro = this->radius_circumscribed_circle();
    if (r > ro)
        return false;

    double ri = 0.5 * sqrt(3.0) * ro;
    if (r <= ri)
        return true;

    // NOTE: Old code used
    //    xl = sqrt(ro^2 - ri^2)
    // where `ro` is the radius of the circumscribing circle and `ri` is
    // the radius of the inscribing circle. But this is equivalent to
    //    xl = 0.5 * ro

    double xl = 0.5 * this->radius_circumscribed_circle();
    double y1, y2;
    if (xl < x && x <= ro)
    {
        y1 = sqrt(3.0) * (x - ro);
        y2 = -y1;
        // if (y1 <= y && y <= y2) return true;
        return (y1 <= y && y <= y2);
    }
    else if (-xl <= x && x <= xl)
    {
        return (-ri <= y && y <= ri);
    }
    else if (-ro <= x && x < -xl)
    {
        y1 = sqrt(3.0) * (x + ro);
        y2 = -y1;
        return (y2 <= y && y <= y1);
    }

    return false;
}

aperture_ptr Hexagon::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<Hexagon>(*this);
}

IrregularTriangle::IrregularTriangle(double x1, double y1,
                                     double x2, double y2,
                                     double x3, double y3)
    : Aperture(ApertureType::IRREGULAR_TRIANGLE),
      x1(x1), y1(y1),
      x2(x2), y2(y2),
      x3(x3), y3(y3)
{
}

double IrregularTriangle::aperture_area() const
{
    double v11 = this->x1 - this->x2;
    double v12 = this->y1 - this->y2;
    double v21 = this->x3 - this->x2;
    double v22 = this->y3 - this->y2;

    double v1m = sqrt(v11 * v11 + v12 * v12);
    double v2m = sqrt(v21 * v21 + v22 * v22);

    double theta = acos((v11 * v21 + v12 * v22) / (v1m * v2m));
    double area = 0.5 * v1m * v2m * sin(theta);

    return area;
}

double IrregularTriangle::diameter_circumscribed_circle() const
{
    // TODO: Not sure this is exact. Is that a problem?
    double xmax = std::max(std::max(x1, x2), x3);
    double ymax = std::max(std::max(y1, y2), y3);
    double xmin = std::min(std::min(x1, x2), x3);
    double ymin = std::min(std::min(y1, y2), y3);
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    return sqrt(dx * dx + dy * dy);
}

bool IrregularTriangle::is_in(double x, double y) const
{
    return intri(x1, y1, x2, y2, x3, y3, x, y);
}

aperture_ptr IrregularTriangle::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<IrregularTriangle>(*this);
}

IrregularQuadrilateral::IrregularQuadrilateral(double x1, double y1,
                                               double x2, double y2,
                                               double x3, double y3,
                                               double x4, double y4)
    : Aperture(ApertureType::IRREGULAR_QUADRILATERAL),
      x1(x1), y1(y1),
      x2(x2), y2(y2),
      x3(x3), y3(y3),
      x4(x4), y4(y4)
{
}

double IrregularQuadrilateral::aperture_area() const
{
    double v11 = this->x1 - this->x2;
    double v12 = this->y1 - this->y2;
    double v21 = this->x3 - this->x2;
    double v22 = this->y3 - this->y2;
    double v31 = this->x3 - this->x4;
    double v32 = this->y3 - this->y4;
    double v41 = this->x1 - this->x4;
    double v42 = this->y1 - this->y4;

    double v1m = sqrt(v11 * v11 + v12 * v12);
    double v2m = sqrt(v21 * v21 + v22 * v22);
    double v3m = sqrt(v31 * v31 + v32 * v32);
    double v4m = sqrt(v41 * v41 + v42 * v42);

    double theta1 = acos((v11 * v21 + v12 * v22) / (v1m * v2m));
    double theta2 = acos((v31 * v41 + v32 * v42) / (v3m * v4m));

    double area = 0.5 * (v1m * v2m * sin(theta1) + v3m * v4m * sin(theta2));
    return area;
}

double IrregularQuadrilateral::diameter_circumscribed_circle() const
{
    // TODO: Not sure this is exact. Is that a problem?
    double xmax = std::max(std::max(x1, x2), std::max(x3, x4));
    double ymax = std::max(std::max(y1, y2), std::max(y3, y4));
    double xmin = std::min(std::min(x1, x2), std::min(x3, x4));
    double ymin = std::min(std::min(y1, y2), std::min(y3, y4));
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    return sqrt(dx * dx + dy * dy);
}

bool IrregularQuadrilateral::is_in(double x, double y) const
{
    return inquad(x1, y1, x2, y2, x3, y3, x4, y4, x, y);
}

aperture_ptr IrregularQuadrilateral::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<IrregularQuadrilateral>(*this);
}

Rectangle::Rectangle(double xlen, double ylen)
    : Aperture(ApertureType::RECTANGLE),
      x_length(xlen),
      y_length(ylen)
{
    // Default to rectangle centered at the origin.
    this->x_coord = -0.5 * this->x_length;
    this->y_coord = -0.5 * this->y_length;
    return;
}

double Rectangle::aperture_area() const
{
    return this->x_length * this->y_length;
}

double Rectangle::diameter_circumscribed_circle() const
{
    return sqrt(x_length * x_length + y_length * y_length);
}

aperture_ptr Rectangle::make_copy() const
{
    // Invokes the implicit copy constructor
    return make_aperture<Rectangle>(*this);
}

bool Rectangle::is_in(double x, double y) const
{
    double xl = this->x_coord;
    double yl = this->y_coord;
    double xu = xl + this->x_length;
    double yu = yl + this->y_length;
    return (xl <= x && x <= xu && yl <= y && y <= yu);
}

Rectangle::Rectangle(double xlen, double ylen, double xl, double yl)
    : Aperture(ApertureType::RECTANGLE),
      x_length(xlen),
      y_length(ylen),
      x_coord(xl),
      y_coord(yl)
{
}

bool intri(double x1, double y1,
           double x2, double y2,
           double x3, double y3,
           double xt, double yt)
{
    double a = (x1 - xt) * (y2 - yt) - (x2 - xt) * (y1 - yt);
    double b = (x2 - xt) * (y3 - yt) - (x3 - xt) * (y2 - yt);
    double c = (x3 - xt) * (y1 - yt) - (x1 - xt) * (y3 - yt);
    return (std::signbit(a) == std::signbit(b) &&
            std::signbit(b) == std::signbit(c));
    // return (sign(a) == sign(b) && sign(b) == sign(c));
}

bool inquad(double x1, double y1,
            double x2, double y2,
            double x3, double y3,
            double x4, double y4,
            double xt, double yt)
{
    return (intri(x1, y1, x2, y2, x3, y3, xt, yt) ||
            intri(x1, y1, x3, y3, x4, y4, xt, yt));
}

} // namespace SolTrace::Data
