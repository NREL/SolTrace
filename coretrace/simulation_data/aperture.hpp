/**
 * @file aperture.hpp
 * @brief Aperture geometry definitions and implementations
 *
 * Defines various aperture shapes (rectangular, circular, annular, etc.) used to
 * limit the active area of optical elements in ray tracing simulations.
 * Apertures determine which rays are accepted or rejected by optical surfaces.
 *
 * @defgroup apertures Aperture Geometries
 * @{
 */

#ifndef SOLTRACE_APERTURE_H
#define SOLTRACE_APERTURE_H

#include <cmath>
#include <memory>
#include <vector>

#include <glm/vec2.hpp>
#include <nlohmann/json.hpp>

#include "constants.hpp"

namespace SolTrace::Data
{

    // TODO: For apertures that do not include the origin, should the
    // "circumscribing" circle be centered at the origin? Or should it
    // be the actual circumscribed circle.

    // Add to ApertureTypeMap when adding to ApertureType
    enum ApertureType
    {
        ANNULUS,
        CIRCLE,
        HEXAGON,
        RECTANGLE,
        EQUILATERAL_TRIANGLE,
        SINGLE_AXIS_CURVATURE_SECTION,
        IRREGULAR_TRIANGLE,
        IRREGULAR_QUADRILATERAL,
        APERTURE_UNKNOWN
    };

    inline const std::map<ApertureType, std::string> ApertureTypeMap =
        {
            {ApertureType::ANNULUS, "ANNULUS"},
            {ApertureType::CIRCLE, "CIRCLE"},
            {ApertureType::HEXAGON, "HEXAGON"},
            {ApertureType::RECTANGLE, "RECTANGLE"},
            {ApertureType::EQUILATERAL_TRIANGLE, "EQUILATERAL_TRIANGLE"},
            {ApertureType::SINGLE_AXIS_CURVATURE_SECTION, "SINGLE_AXIS_CURVATURE_SECTION"},
            {ApertureType::IRREGULAR_TRIANGLE, "IRREGULAR_TRIANGLE"},
            {ApertureType::IRREGULAR_QUADRILATERAL, "IRREGULAR_QUADRILATERAL"},
            {ApertureType::APERTURE_UNKNOWN, "APERTURE_UNKNOWN"}};

    struct Aperture;
    using aperture_ptr = std::shared_ptr<Aperture>;

    /**
     * @brief Factory function to create aperture objects
     * @tparam A Aperture type to create
     * @tparam Args Constructor argument types
     * @param args Constructor arguments
     * @return Shared pointer to the created aperture
     */
    template <typename A, typename... Args>
    inline auto make_aperture(Args &&...args)
    {
        return std::make_shared<A>(std::forward<Args>(args)...);
    }

    struct Aperture
    {
    public:
        ApertureType my_type;

        /**
         * @brief Constructor for base aperture
         * @param type The aperture type enumeration
         */
        Aperture(ApertureType type) : my_type(type) {}

        virtual ~Aperture() {}

        /**
         * @brief Factory method to create apertures from type and parameters
         * @param type The aperture type to create
         * @param args Vector of parameters for the aperture
         * @return Shared pointer to the created aperture
         */
        static aperture_ptr make_aperture_from_type(ApertureType type,
                                                    const std::vector<double> &args);

        /**
         * @brief Factory method to create apertures from json
         * @param jnode the json containing necessary parameters for each aperture type
         * @return Shared pointer to the created aperture
         */
        static aperture_ptr make_aperture_from_json(const nlohmann::ordered_json &jnode);

        /**
         * @brief Get the aperture type
         * @return The aperture type enumeration
         */
        inline ApertureType get_type() const
        {
            return my_type;
        }

        /**
         * @brief Get radius of circumscribed circle
         * @return Radius of the smallest circle that contains the aperture
         */
        virtual inline double radius_circumscribed_circle() const
        {
            return 0.5 * this->diameter_circumscribed_circle();
        }

        /**
         * @brief Calculate the aperture area
         * @return Area of the aperture in square units
         */
        virtual double aperture_area() const = 0;

        /**
         * @brief Get diameter of circumscribed circle
         * @return Diameter of the smallest circle that contains the aperture
         */
        virtual double diameter_circumscribed_circle() const = 0;

        /**
         * @brief Triangulate the aperture shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const = 0;

        /**
         * @brief Test if a point is inside the aperture
         * @param x X coordinate of the test point
         * @param y Y coordinate of the test point
         * @return True if point is inside aperture, false otherwise
         */
        virtual bool is_in(double x, double y) const = 0;

        /**
         * @brief Create a copy of this aperture
         * @return Shared pointer to a copy of this aperture
         */
        virtual aperture_ptr make_copy() const = 0;

        virtual void bounding_box(double &xmin, double &xmax,
                                  double &ymin, double &ymax) const = 0;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const = 0;

        /**
         * @brief Get the aperture type string
         * @return The aperture type string
         */
        inline std::string get_type_string() const
        {
            switch (my_type)
            {
            case ANNULUS:
                return "Annulus";
            case CIRCLE:
                return "Circle";
            case HEXAGON:
                return "Hexagon";
            case RECTANGLE:
                return "Rectangle";
            case EQUILATERAL_TRIANGLE:
                return "Regular Triangle";
            case SINGLE_AXIS_CURVATURE_SECTION:
                return "Single Axis Curvature";
            case IRREGULAR_TRIANGLE:
                return "Triangle";
            case IRREGULAR_QUADRILATERAL:
                return "Quad";
            case APERTURE_UNKNOWN:
                return "Unknown";
            }
            return "Unknown";
        }

    protected:
        struct Point
        {
        public:
            double x;
            double y;
            Point(double ix, double iy) : x(ix), y(iy) {}
            bool operator==(const Point &p) const
            {
                return x == p.x && y == p.y;
            }
        };
        struct Triangle
        {
        public:
            Point a;
            Point b;
            Point c;
            Triangle(Point ia, Point ib, Point ic) : a(ia), b(ib), c(ic) {}
        };
        /**
         * @brief Compute the mipoint between to points
         * @return The midpoint
         */
        Point midpoint(const Point &v0, const Point &v1) const;
        /**
         * @brief Recursively Subdivide a triangle by midpoints
         * @param tri The triangle to subdivide
         * @param n number of subdivisions
         * @return vector of triangles
         */
        std::vector<Triangle> subdivide(Triangle tri, int n) const;
        /**
         * @brief Return the index of the point in the vector, adding if it
         * doesn't exist
         * @param v The vector of points
         * @param p The point of interest
         * @return index of p in v
         */
        int index_of(std::vector<Point> &v, const Point &p) const;
        /**
         * @brief Convert a list of Triangles in indexed (flattened) faceset
         * @param triangles The list to convert
         * @return indexed faceset
         */
        std::tuple<std::vector<double>, std::vector<int>>
        indexed_triangles(const std::vector<Triangle> &triangles) const;
    };

    struct Annulus : public Aperture
    {
        double inner_radius;
        double outer_radius;
        double arc_angle;

        // Annulus()
        //     : Aperture(ANNULUS),
        //       inner_radius(0.0), outer_radius(0.0), arc_angle(0.0)
        // {
        // }

        /**
         * @brief Constructor for annulus aperture
         * @param ri Inner radius of the annulus
         * @param ro Outer radius of the annulus
         * @param arc Arc angle in radians (2*pi for full annulus)
         */
        Annulus(double ri, double ro, double arc)
            : Aperture(ApertureType::ANNULUS),
              inner_radius(ri),
              outer_radius(ro),
              arc_angle(arc)
        {
        }

        /**
         * @brief Json-based constructor for annulus aperture
         * @param jnode contains ri, ro, and arc angle
         */
        Annulus(const nlohmann::ordered_json &jnode);

        virtual ~Annulus() {}

        /**
         * @brief Calculate annulus aperture area
         * @return Area of the annular aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for annulus
         * @return Diameter of outer circle
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside annulus aperture
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within annular region
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of annulus aperture
         * @return Shared pointer to annulus copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the annulus shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    struct Circle : public Aperture
    {
        double diameter;

        /**
         * @brief Constructor for circular aperture
         * @param d Diameter of the circle
         */
        Circle(double d) : Aperture(ApertureType::CIRCLE), diameter(d) {}

        /**
         * @brief Json-based constructor for circular aperture
         * @param jnode contains diameter
         */
        Circle(const nlohmann::ordered_json &jnode);

        virtual ~Circle() {}

        /**
         * @brief Calculate circular aperture area
         * @return Area of the circular aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for circle
         * @return Circle diameter (same as aperture diameter)
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside circular aperture
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within circle
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of circular aperture
         * @return Shared pointer to circle copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the circle shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    struct EquilateralTriangle : public Aperture
    {
        double circumscribe_diameter;
        // EquilateralTriangle() : Aperture(EQUILATERAL_TRIANGLE),
        //                         circumscribe_diameter(0.0)
        // {
        // }

        /**
         * @brief Constructor for equilateral triangle aperture
         * @param cd Diameter of circumscribed circle
         */
        EquilateralTriangle(double cd)
            : Aperture(ApertureType::EQUILATERAL_TRIANGLE),
              circumscribe_diameter(cd) { }

        /**
         * @brief Json-based constructor for equilateral triangle aperture
         * @param jnode contains diameter of circumscribed circle
         */
        EquilateralTriangle(const nlohmann::ordered_json &jnode);

        virtual ~EquilateralTriangle() { }

        /**
         * @brief Calculate equilateral triangle aperture area
         * @return Area of the triangular aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for triangle
         * @return Diameter of circumscribed circle
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside triangular aperture
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within triangle
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of triangular aperture
         * @return Shared pointer to triangle copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the triangle shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    struct Hexagon : public Aperture
    {
        double circumscribe_diameter;

        // Hexagon() : Aperture(HEXAGON), circumscribe_diameter(0.0) {}

        /**
         * @brief Constructor for hexagonal aperture
         * @param d Diameter of circumscribed circle
         */
        Hexagon(double d)
            : Aperture(ApertureType::HEXAGON),
              circumscribe_diameter(d) {}

        /**
         * @brief Json-based constructor for hexagonal aperture
         * @param jnode contains diameter of circumscribed circle
         */
        Hexagon(const nlohmann::ordered_json &jnode);

        virtual ~Hexagon() {}

        /**
         * @brief Calculate hexagonal aperture area
         * @return Area of the hexagonal aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for hexagon
         * @return Diameter of circumscribed circle
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside hexagonal aperture
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within hexagon
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of hexagonal aperture
         * @return Shared pointer to hexagon copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the hexagon shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    class Rectangle : public Aperture
    {
        glm::dvec2 m_length;
        // NOTE: The point (x_coord, y_coord) gives the location of the
        // lower left hand corner of the rectangle in the xy-plane.
        glm::dvec2 m_coord;

        glm::dvec2 m_cached_range;

        void update_cached() { m_cached_range = m_coord + m_length; }

    public:
        /**
         * @brief Constructor for centered rectangular aperture
         * @param xlen Length in x direction
         * @param ylen Length in y direction
         */
        Rectangle(double xlen, double ylen);

        /**
         * @brief Constructor for positioned rectangular aperture
         * @param xlen Length in x direction
         * @param ylen Length in y direction
         * @param xl X coordinate of lower-left corner
         * @param yl Y coordinate of lower-left corner
         */
        Rectangle(double xlen, double ylen, double xl, double yl);

        /**
         * @brief Json-based constructor for rectangular aperture
         * @param jnode contains xlen, ylen, xl, yl
         */
        Rectangle(const nlohmann::ordered_json &jnode);

        virtual ~Rectangle() {}

        double x_length() const { return m_length.x; }
        double y_length() const { return m_length.y; }
        double x_coord() const { return m_coord.x; }
        double y_coord() const { return m_coord.y; }

        void set_x_length(double x_length)
        {
            m_length.x = x_length;
            update_cached();
        }
        void set_y_length(double y_length)
        {
            m_length.y = y_length;
            update_cached();
        }
        void set_x_coord(double x_coord)
        {
            m_coord.x = x_coord;
            update_cached();
        }
        void set_y_coord(double y_coord)
        {
            m_coord.y = y_coord;
            update_cached();
        }

        /**
         * @brief Calculate rectangular aperture area
         * @return Area of the rectangular aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for rectangle
         * @return Diagonal length of rectangle
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside rectangular aperture
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within rectangle bounds
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of rectangular aperture
         * @return Shared pointer to rectangle copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the rectangle shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    // struct SingleAxisCurvatureSection : public Aperture
    // {
    //     // TODO: Implement this?
    // };

    struct IrregularTriangle : public Aperture
    {
        // Locations of the 3 vertices
        double x1;
        double y1;
        double x2;
        double y2;
        double x3;
        double y3;

        /**
         * @brief Constructor for irregular triangle aperture
         * @param x1 X coordinate of vertex 1
         * @param y1 Y coordinate of vertex 1
         * @param x2 X coordinate of vertex 2
         * @param y2 Y coordinate of vertex 2
         * @param x3 X coordinate of vertex 3
         * @param y3 Y coordinate of vertex 3
         */
        IrregularTriangle(double x1, double y1,
                          double x2, double y2,
                          double x3, double y3);

        /**
         * @brief Json-based constructor for irregular triangle aperture
         * @param jnode contains x1, y1, x2, y2, x3, y3
         */
        IrregularTriangle(const nlohmann::ordered_json &jnode);

        ~IrregularTriangle() {}

        /**
         * @brief Calculate irregular triangle aperture area
         * @return Area of the triangular aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for triangle
         * @return Diameter of smallest circle containing triangle
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside irregular triangle
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within triangle
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of irregular triangle aperture
         * @return Shared pointer to triangle copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the triangle shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    struct IrregularQuadrilateral : public Aperture
    {
        // Locations of the 4 vertices
        double x1;
        double y1;
        double x2;
        double y2;
        double x3;
        double y3;
        double x4;
        double y4;

        /**
         * @brief Constructor for irregular quadrilateral aperture
         * @param x1 X coordinate of vertex 1
         * @param y1 Y coordinate of vertex 1
         * @param x2 X coordinate of vertex 2
         * @param y2 Y coordinate of vertex 2
         * @param x3 X coordinate of vertex 3
         * @param y3 Y coordinate of vertex 3
         * @param x4 X coordinate of vertex 4
         * @param y4 Y coordinate of vertex 4
         */
        IrregularQuadrilateral(double x1, double y1,
                               double x2, double y2,
                               double x3, double y3,
                               double x4, double y4);

        /**
         * @brief Json-based constructor for irregular quadrilateral aperture
         * @param jnode contains x1, y1, x2, y2, x3, y3, x4, y4
         */
        IrregularQuadrilateral(const nlohmann::ordered_json &jnode);

        ~IrregularQuadrilateral() {}

        /**
         * @brief Calculate irregular quadrilateral aperture area
         * @return Area of the quadrilateral aperture
         */
        virtual double aperture_area() const override;

        /**
         * @brief Get diameter of circumscribed circle for quadrilateral
         * @return Diameter of smallest circle containing quadrilateral
         */
        virtual double diameter_circumscribed_circle() const override;

        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override;

        /**
         * @brief Test if point is inside irregular quadrilateral
         * @param x X coordinate
         * @param y Y coordinate
         * @return True if point is within quadrilateral
         */
        virtual bool is_in(double x, double y) const override;

        /**
         * @brief Create copy of irregular quadrilateral aperture
         * @return Shared pointer to quadrilateral copy
         */
        virtual aperture_ptr make_copy() const override;

        /**
         * @brief Write aperture parameters to json
         * @param jnode JSON node
         */
        virtual void write_json(nlohmann::ordered_json &jnode) const override;

        /**
         * @brief Triangulate the quad shape
         * @return Tuple of 2D vertices and triangle indices
         */
        virtual std::tuple<std::vector<double>, std::vector<int>>
        triangulation() const override;
    };

    /**
     * @brief Test if point is inside triangle defined by three vertices
     * @param x1 X coordinate of vertex 1
     * @param y1 Y coordinate of vertex 1
     * @param x2 X coordinate of vertex 2
     * @param y2 Y coordinate of vertex 2
     * @param x3 X coordinate of vertex 3
     * @param y3 Y coordinate of vertex 3
     * @param xt X coordinate of test point
     * @param yt Y coordinate of test point
     * @return True if point is inside triangle
     */
    bool intri(double x1, double y1,
               double x2, double y2,
               double x3, double y3,
               double xt, double yt);

    /**
     * @brief Test if point is inside quadrilateral defined by four vertices
     * @param x1 X coordinate of vertex 1
     * @param y1 Y coordinate of vertex 1
     * @param x2 X coordinate of vertex 2
     * @param y2 Y coordinate of vertex 2
     * @param x3 X coordinate of vertex 3
     * @param y3 Y coordinate of vertex 3
     * @param x4 X coordinate of vertex 4
     * @param y4 Y coordinate of vertex 4
     * @param xt X coordinate of test point
     * @param yt Y coordinate of test point
     * @return True if point is inside quadrilateral
     */
    bool inquad(double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                double x4, double y4,
                double xt, double yt);

} // namespace SolTrace::Data

/**
 * @}
 */

#endif
