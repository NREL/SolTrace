
#include "aperture.hpp"
#include "simdata_io.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

// #include <iostream>

#include "constants.hpp"

namespace SolTrace::Data
{

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

    aperture_ptr Aperture::make_aperture_from_json(const nlohmann::ordered_json &jnode)
    {
        if (!jnode.contains("aperture_type"))
            throw std::invalid_argument("Missing aperture_type");
        std::string type_str = jnode.at("aperture_type");
        ApertureType aperture_type = get_enum_from_string(type_str, ApertureTypeMap, ApertureType::APERTURE_UNKNOWN);
        switch (aperture_type)
        {
        case ApertureType::ANNULUS:
            return make_aperture<Annulus>(jnode);
        case ApertureType::CIRCLE:
            return make_aperture<Circle>(jnode);
        case ApertureType::HEXAGON:
            return make_aperture<Hexagon>(jnode);
        case ApertureType::RECTANGLE:
            return make_aperture<Rectangle>(jnode);
        case ApertureType::EQUILATERAL_TRIANGLE:
            return make_aperture<EqualateralTriangle>(jnode);
        case ApertureType::IRREGULAR_TRIANGLE:
            return make_aperture<IrregularTriangle>(jnode);
        case ApertureType::IRREGULAR_QUADRILATERAL:
            return make_aperture<IrregularQuadrilateral>(jnode);
        default:
            throw std::invalid_argument("Unsupported aperture_type: " + type_str);
        }
    }

    Aperture::Point Aperture::midpoint(const Aperture::Point &v0,
                                       const Aperture::Point &v1) const
    {
        return Aperture::Point((v0.x + v1.x) / 2, (v0.y + v1.y) / 2);
    }

    std::vector<Aperture::Triangle> Aperture::subdivide(Aperture::Triangle tri,
                                                        int n) const
    {
        Point v0 = tri.a;
        Point v1 = tri.b;
        Point v2 = tri.c;
        Point m01 = midpoint(v0, v1);
        Point m12 = midpoint(v1, v2);
        Point m20 = midpoint(v2, v0);
        std::vector<Triangle> result;
        if (n - 1 == 0)
        {
            result.push_back(Triangle(v0, m01, m20));
            result.push_back(Triangle(m01, v1, m12));
            result.push_back(Triangle(m12, m01, m20));
            result.push_back(Triangle(m12, m20, v2));
            return result;
        }
        auto t1 = subdivide(Triangle(v0, m01, m20), n - 1);
        auto t2 = subdivide(Triangle(v0, m01, m20), n - 1);
        auto t3 = subdivide(Triangle(v0, m01, m20), n - 1);
        auto t4 = subdivide(Triangle(v0, m01, m20), n - 1);
        result.insert(result.end(), t1.begin(), t1.end());
        result.insert(result.end(), t2.begin(), t2.end());
        result.insert(result.end(), t3.begin(), t3.end());
        result.insert(result.end(), t4.begin(), t4.end());
        return result;
    }

    int Aperture::index_of(std::vector<Aperture::Point> &v,
                           const Aperture::Point &p) const
    {
        auto it = find(v.begin(), v.end(), p);
        if (it == v.end())
        {
            v.push_back(p);
            return v.size() - 1;
        }
        return std::distance(v.begin(), it);
    }

    std::tuple<std::vector<double>, std::vector<int>> Aperture::indexed_triangles(
        const std::vector<Aperture::Triangle> &triangles) const
    {
        std::vector<int> indices;
        std::vector<Point> points;
        std::vector<double> flattened;
        for (const Triangle &tri : triangles)
        {
            indices.push_back(index_of(points, tri.a));
            indices.push_back(index_of(points, tri.b));
            indices.push_back(index_of(points, tri.c));
        }
        for (const Point &p : points)
        {
            flattened.push_back(p.x);
            flattened.push_back(p.y);
        }
        return std::make_pair(flattened, indices);
    }

    Annulus::Annulus(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::ANNULUS)
    {
        this->inner_radius = jnode.at("inner_radius");
        this->outer_radius = jnode.at("outer_radius");
        this->arc_angle = jnode.at("arc_angle");
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

    void Annulus::bounding_box(double &xmin,
                               double &xmax,
                               double &ymin,
                               double &ymax) const
    {
        xmin = -this->outer_radius;
        xmax = this->outer_radius;
        ymin = -this->outer_radius;
        ymax = this->outer_radius;
        return;
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

    std::tuple<std::vector<double>, std::vector<int>>
    Annulus::triangulation() const
    {
        const int resolution = 32;
        std::vector<double> verts;
        std::vector<int> indices;
        for (int i = 0; i <= resolution; i++)
        {
            const double u = i / resolution * PI * 2;
            verts.push_back(inner_radius * std::cos(u));
            verts.push_back(inner_radius * std::sin(u));
            verts.push_back(outer_radius * std::cos(u));
            verts.push_back(outer_radius * std::sin(u));
        }
        for (int i = 0; i < resolution - 3; i += 2)
        {
            const int a = i;
            const int b = i + 1;
            const int c = i + 2;
            const int d = i + 3;
            // Generate two triangles for each quad in the mesh
            // Adjust order to be counter-clockwise
            indices.push_back(a);
            indices.push_back(d);
            indices.push_back(b);
            indices.push_back(b);
            indices.push_back(d);
            indices.push_back(c);
        }
        return std::make_tuple(verts, indices);
    }

    aperture_ptr Annulus::make_copy() const
    {
        // Invokes the implicit copy constructor
        return make_aperture<Annulus>(*this);
    }

    void Annulus::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::ANNULUS;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["inner_radius"] = this->inner_radius;
        jnode["outer_radius"] = this->outer_radius;
        jnode["arc_angle"] = this->arc_angle;
    }

    Circle::Circle(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::CIRCLE)
    {
        this->diameter = jnode.at("diameter");
    }

    double Circle::aperture_area() const
    {
        return 0.25 * PI * this->diameter * this->diameter;
    }

    double Circle::diameter_circumscribed_circle() const
    {
        return this->diameter;
    }

    void Circle::bounding_box(double &xmin,
                              double &xmax,
                              double &ymin,
                              double &ymax) const
    {
        double r = this->radius_circumscribed_circle();
        xmin = -r;
        xmax = r;
        ymin = -r;
        ymax = r;
        return;
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

    void Circle::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::CIRCLE;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["diameter"] = this->diameter;
    }

    std::tuple<std::vector<double>, std::vector<int>>
    Circle::triangulation() const
    {
        // Using a fixed Delaunay triangulation of the unit circle
        std::vector<double> verts = {
            1.00000000e+00, 0.00000000e+00, 9.51056516e-01, 3.09016994e-01,
            8.09016994e-01, 5.87785252e-01, 5.87785252e-01, 8.09016994e-01,
            3.09016994e-01, 9.51056516e-01, 6.12323400e-17, 1.00000000e+00,
            -3.09016994e-01, 9.51056516e-01, -5.87785252e-01, 8.09016994e-01,
            -8.09016994e-01, 5.87785252e-01, -9.51056516e-01, 3.09016994e-01,
            -1.00000000e+00, 1.22464680e-16, -9.51056516e-01, -3.09016994e-01,
            -8.09016994e-01, -5.87785252e-01, -5.87785252e-01, -8.09016994e-01,
            -3.09016994e-01, -9.51056516e-01, -1.83697020e-16, -1.00000000e+00,
            3.09016994e-01, -9.51056516e-01, 5.87785252e-01, -8.09016994e-01,
            8.09016994e-01, -5.87785252e-01, 9.51056516e-01, -3.09016994e-01,
            0.00000000e+00, 0.00000000e+00, -5.00000000e-01, -5.00000000e-01,
            -5.00000000e-01, 0.00000000e+00, -5.00000000e-01, 5.00000000e-01,
            0.00000000e+00, -5.00000000e-01, 0.00000000e+00, 5.00000000e-01,
            5.00000000e-01, -5.00000000e-01, 5.00000000e-01, 0.00000000e+00,
            5.00000000e-01, 5.00000000e-01};
        std::vector<int> indices = {
            22, 11, 21, 11, 22, 10, 17, 18, 26, 14, 24, 21, 13, 14, 21, 24, 14, 15,
            25, 6, 23, 6, 25, 5, 9, 22, 23, 8, 9, 23, 22, 9, 10, 27, 1, 28,
            1, 27, 0, 19, 27, 26, 18, 19, 26, 27, 19, 0, 4, 25, 28, 3, 4, 28,
            25, 4, 5, 12, 13, 21, 11, 12, 21, 24, 16, 26, 16, 17, 26, 16, 24, 15,
            7, 8, 23, 6, 7, 23, 2, 3, 28, 1, 2, 28, 24, 22, 21, 22, 24, 20,
            24, 27, 20, 27, 24, 26, 25, 22, 20, 22, 25, 23, 27, 25, 20, 25, 27, 28};
        // scale from the unit cirle to our circle
        std::transform(
            verts.begin(), verts.end(), verts.begin(), [this](double element)
            { return element *= this->diameter / 2.0; });
        return std::make_tuple(verts, indices);
    }

    EqualateralTriangle::EqualateralTriangle(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::EQUILATERAL_TRIANGLE)
    {
        this->circumscribe_diameter = jnode.at("circumscribe_diameter");
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

    void EqualateralTriangle::bounding_box(double &xmin,
                                           double &xmax,
                                           double &ymin,
                                           double &ymax) const
    {
        double r = this->radius_circumscribed_circle();
        double side = sqrt(3.0) * r;
        double h = 1.5 * r;
        xmin = -0.5 * side;
        xmax = 0.5 * side;
        ymin = r - h;
        ymax = r;
        return;
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

    std::tuple<std::vector<double>, std::vector<int>>
    EqualateralTriangle::triangulation() const
    {
        double r = circumscribe_diameter / 2.0;
        Triangle tri(Point(0, r),
                     Point(r * cos(-PI / 6.0), r * sin(-PI / 6.0)),
                     Point(r * cos(7 * PI / 6.0), r * sin(7 * PI / 6.0)));
        return indexed_triangles(subdivide(tri, 3));
    }

    aperture_ptr EqualateralTriangle::make_copy() const
    {
        // Invokes the implicit copy constructor
        return make_aperture<EqualateralTriangle>(*this);
    }

    void EqualateralTriangle::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::EQUILATERAL_TRIANGLE;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["circumscribe_diameter"] = this->circumscribe_diameter;
    }

    Hexagon::Hexagon(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::HEXAGON)
    {
        this->circumscribe_diameter = jnode.at("circumscribe_diameter");
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

    void Hexagon::bounding_box(double &xmin,
                               double &xmax,
                               double &ymin,
                               double &ymax) const
    {
        double r = this->radius_circumscribed_circle();
        double apothem = 0.5 * r * sqrt(3.0);
        xmin = -r;
        xmax = r;
        ymin = -apothem;
        ymax = apothem;
        return;
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

    std::tuple<std::vector<double>, std::vector<int>>
    Hexagon::triangulation() const
    {
        double r = circumscribe_diameter / 2.0;
        std::vector<Triangle> t0 =
            subdivide(Triangle(Point(0, 0),
                               Point(r * cos(PI / 3.0), r * sin(PI / 3.0)),
                               Point(r, 0)),
                      2);
        std::vector<Triangle> t1 = subdivide(
            Triangle(Point(0, 0),
                     Point(r * cos(2 * PI / 3.0), r * sin(2 * PI / 3.0)),
                     Point(r * cos(PI / 3.0), r * sin(PI / 3.0))),
            2);
        std::vector<Triangle> t2 = subdivide(
            Triangle(Point(0, 0),
                     Point(-r, 0),
                     Point(r * cos(2 * PI / 3.0), r * sin(2 * PI / 3.0))),
            2);
        std::vector<Triangle> t3 = subdivide(
            Triangle(Point(0, 0),
                     Point(-r, 0),
                     Point(r * cos(4 * PI / 3.0), r * sin(4 * PI / 3.0))),
            2);
        std::vector<Triangle> t4 = subdivide(
            Triangle(Point(0, 0),
                     Point(r * cos(5 * PI / 3.0), r * sin(5 * PI / 3.0)),
                     Point(r * cos(4 * PI / 3.0), r * sin(4 * PI / 3.0))),
            2);
        std::vector<Triangle> t5 = subdivide(
            Triangle(Point(0, 0),
                     Point(r, 0),
                     Point(r * cos(5 * PI / 3.0), r * sin(5 * PI / 3.0))),
            2);
        std::vector<Triangle> triangles;
        triangles.insert(triangles.end(), t0.begin(), t0.end());
        triangles.insert(triangles.end(), t1.begin(), t1.end());
        triangles.insert(triangles.end(), t2.begin(), t2.end());
        triangles.insert(triangles.end(), t3.begin(), t3.end());
        triangles.insert(triangles.end(), t4.begin(), t4.end());
        triangles.insert(triangles.end(), t5.begin(), t5.end());
        return indexed_triangles(triangles);
    }

    aperture_ptr Hexagon::make_copy() const
    {
        // Invokes the implicit copy constructor
        return make_aperture<Hexagon>(*this);
    }

    void Hexagon::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::HEXAGON;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["circumscribe_diameter"] = this->circumscribe_diameter;
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

    IrregularTriangle::IrregularTriangle(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::IRREGULAR_TRIANGLE)
    {
        this->x1 = jnode.at("x1");
        this->y1 = jnode.at("y1");
        this->x2 = jnode.at("x2");
        this->y2 = jnode.at("y2");
        this->x3 = jnode.at("x3");
        this->y3 = jnode.at("y3");
    }

    std::tuple<std::vector<double>, std::vector<int>>
    IrregularTriangle::triangulation() const
    {
        Triangle tri(Point(x1, y1), Point(x2, y2), Point(x3, y3));
        return indexed_triangles(subdivide(tri, 3));
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

    void IrregularTriangle::bounding_box(double &xmin,
                                         double &xmax,
                                         double &ymin,
                                         double &ymax) const
    {
        xmin = std::min(x1, std::min(x2, x3));
        xmax = std::max(x1, std::max(x2, x3));
        ymin = std::min(y1, std::min(y2, y3));
        ymax = std::max(y1, std::max(y2, y3));
        return;
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

    void IrregularTriangle::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::IRREGULAR_TRIANGLE;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["x1"] = this->x1;
        jnode["y1"] = this->y1;
        jnode["x2"] = this->x2;
        jnode["y2"] = this->y2;
        jnode["x3"] = this->x3;
        jnode["y3"] = this->y3;
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

    IrregularQuadrilateral::IrregularQuadrilateral(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::IRREGULAR_QUADRILATERAL)
    {
        this->x1 = jnode.at("x1");
        this->y1 = jnode.at("y1");
        this->x2 = jnode.at("x2");
        this->y2 = jnode.at("y2");
        this->x3 = jnode.at("x3");
        this->y3 = jnode.at("y3");
        this->x4 = jnode.at("x4");
        this->y4 = jnode.at("y4");
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

    std::tuple<std::vector<double>, std::vector<int>>
    IrregularQuadrilateral::triangulation() const
    {
        std::vector<Triangle> t0 =
            subdivide(Triangle(Point(x1, y1), Point(x3, y3), Point(x2, y2)), 2);
        std::vector<Triangle> t1 =
            subdivide(Triangle(Point(x1, y1), Point(x4, y4), Point(x3, y3)), 2);
        std::vector<Triangle> triangles;
        triangles.insert(triangles.end(), t0.begin(), t0.end());
        triangles.insert(triangles.end(), t1.begin(), t1.end());
        return indexed_triangles(triangles);
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

    void IrregularQuadrilateral::bounding_box(double &xmin,
                                              double &xmax,
                                              double &ymin,
                                              double &ymax) const
    {
        xmin = std::min(std::min(x1, x2), std::min(x3, x4));
        xmax = std::max(std::max(x1, x2), std::max(x3, x4));
        ymin = std::min(std::min(y1, y2), std::min(y3, y4));
        ymax = std::max(std::max(y1, y2), std::max(y3, y4));
        return;
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

    void IrregularQuadrilateral::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::IRREGULAR_QUADRILATERAL;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["x1"] = this->x1;
        jnode["y1"] = this->y1;
        jnode["x2"] = this->x2;
        jnode["y2"] = this->y2;
        jnode["x3"] = this->x3;
        jnode["y3"] = this->y3;
        jnode["x4"] = this->x4;
        jnode["y4"] = this->y4;
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

    Rectangle::Rectangle(double xlen, double ylen, double xl, double yl)
        : Aperture(ApertureType::RECTANGLE),
          x_length(xlen),
          y_length(ylen),
          x_coord(xl),
          y_coord(yl)
    {
    }

    Rectangle::Rectangle(const nlohmann::ordered_json &jnode)
        : Aperture(ApertureType::RECTANGLE)
    {
        this->x_length = jnode.at("x_length");
        this->y_length = jnode.at("y_length");
        this->x_coord = jnode.at("x_coord");
        this->y_coord = jnode.at("y_coord");
    }

    double Rectangle::aperture_area() const
    {
        return this->x_length * this->y_length;
    }

    double Rectangle::diameter_circumscribed_circle() const
    {
        return sqrt(x_length * x_length + y_length * y_length);
    }

    void Rectangle::bounding_box(double &xmin,
                                 double &xmax,
                                 double &ymin,
                                 double &ymax) const
    {
        xmin = this->x_coord;
        xmax = this->x_coord + this->x_length;
        ymin = this->y_coord;
        ymax = this->y_coord + this->y_length;
        return;
    }

    bool Rectangle::is_in(double x, double y) const
    {
        double xl = this->x_coord;
        double yl = this->y_coord;
        double xu = xl + this->x_length;
        double yu = yl + this->y_length;
        return (xl <= x && x <= xu && yl <= y && y <= yu);
    }

    aperture_ptr Rectangle::make_copy() const
    {
        // Invokes the implicit copy constructor
        return make_aperture<Rectangle>(*this);
    }

    void Rectangle::write_json(nlohmann::ordered_json &jnode) const
    {
        ApertureType type = ApertureType::RECTANGLE;
        jnode["aperture_type"] = ApertureTypeMap.at(type);
        jnode["x_length"] = this->x_length;
        jnode["y_length"] = this->y_length;
        jnode["x_coord"] = this->x_coord;
        jnode["y_coord"] = this->y_coord;
    }

    std::tuple<std::vector<double>, std::vector<int>>
    Rectangle::triangulation() const
    {
        const int segments = 5;
        std::vector<double> verts;
        std::vector<int> indices;
        for (int i = 0; i <= segments; ++i)
        {
            for (int j = 0; j <= segments; ++j)
            {
                const double x = i * x_length / segments + x_coord;
                const double y = j * y_length / segments + y_coord;
                verts.push_back(x);
                verts.push_back(y);
            }
        }
        for (int i = 0; i < segments; ++i)
        {
            for (int j = 0; j < segments; ++j)
            {
                const int a = (segments + 1) * i + j;
                const int c = (segments + 1) * (i + 1) + j;
                const int d = (segments + 1) * (i + 1) + j + 1;
                const int b = (segments + 1) * i + j + 1;
                // Generate two triangles for each quad in the mesh
                // Adjust order to be counter-clockwise
                indices.push_back(a);
                indices.push_back(c);
                indices.push_back(b);
                indices.push_back(b);
                indices.push_back(c);
                indices.push_back(d);
            }
        }
        return make_pair(verts, indices);
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
