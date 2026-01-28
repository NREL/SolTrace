#include "bbox_calculator.hpp"

#include <algorithm>

#include <matvec.hpp>
#include <native_runner_types.hpp>

namespace SolTrace::EmbreeRunner
{
    using SolTrace::Data::MatrixVectorMult_generic;

    using SolTrace::NativeRunner::TElement;
    using SolTrace::NativeRunner::telement_ptr;
    using SolTrace::NativeRunner::TStage;
    using SolTrace::NativeRunner::tstage_ptr;

    // float get_absolute_minmax(const float values[], int size, bool is_max)
    // {
    //     float minmax_abs = std::abs(values[0]);

    //     if (is_max)
    //     {
    //         for (int i = 1; i < size; i++)
    //         {
    //             minmax_abs = std::max(std::abs(values[i]),
    //                                   minmax_abs);
    //         }
    //     }
    //     else
    //     {
    //         for (int i = 1; i < size; i++)
    //         {
    //             minmax_abs = std::min(std::abs(values[i]),
    //                                   minmax_abs);
    //         }
    //     }

    //     return minmax_abs;
    // }

    // void process_zernike_bounds(TElement* st_element, float x_minmax[2], float y_minmax[2],
    // 	float(&z_minmax)[2])
    // {
    // 	double z, x, y;

    // 	float& z_min = z_minmax[0];
    // 	float& z_max = z_minmax[1];

    // 	z_min = std::numeric_limits<double>::infinity();
    // 	z_max = -std::numeric_limits<double>::infinity();

    // 	float x_range[3] = { x_minmax[0], x_minmax[1], 0.f };
    // 	float y_range[3] = { y_minmax[0], y_minmax[1], 0.f };

    // 	for (int xi = 0; xi < 3; ++xi)
    // 	{
    // 		for (int yi = 0; yi < 3; ++yi)
    // 		{
    // 			x = x_range[xi];
    // 			y = y_range[yi];
    // 			EvalMono(x, y, st_element->BCoefficients, st_element->FitOrder, 0.0, 0.0, &z);

    // 			if (z < z_min) z_min = z;
    // 			if (z > z_max) z_max = z;
    // 		}
    // 	}
    // }

    // void process_poly_bounds(TElement* st_element, float x_minmax[2], float y_minmax[2],
    // 	float(&z_minmax)[2])
    // {
    // 	double z, x, y;

    // 	float& z_min = z_minmax[0];
    // 	float& z_max = z_minmax[1];

    // 	z_min = std::numeric_limits<double>::infinity();
    // 	z_max = -std::numeric_limits<double>::infinity();

    // 	float x_range[3] = { x_minmax[0], x_minmax[1], 0.f };
    // 	float y_range[3] = { y_minmax[0], y_minmax[1], 0.f };

    // 	for (int xi = 0; xi < 3; ++xi)
    // 	{
    // 		for (int yi = 0; yi < 3; ++yi)
    // 		{
    // 			x = x_range[xi];
    // 			y = y_range[yi];
    // 			EvalPoly(x, y, st_element->PolyCoeffs, st_element->FitOrder, &z);

    // 			if (z < z_min) z_min = z;
    // 			if (z > z_max) z_max = z;
    // 		}
    // 	}
    // }

    // bool find_spline_extrema(std::vector<double>& xa,
    // 	std::vector<double>& ya,
    // 	std::vector<double>& y2a,
    // 	double xMin, double xMax,
    // 	double& yMin, double& yMax)
    // {
    // 	if (xa.size() < 2 || xa.size() != ya.size() || xa.size() != y2a.size())
    // 		return false;

    // 	int n = xa.size();
    // 	yMin = INFINITY;
    // 	yMax = -INFINITY;

    // 	// Check all interval critical points
    // 	for (int i = 0; i < n - 1; i++) {
    // 		double h = xa[i + 1] - xa[i];
    // 		if (h == 0) continue;

    // 		// The spline derivative equal to zero gives us a quadratic equation
    // 		// Coefficients derived from the derivative expression in splint function
    // 		double A = h * (y2a[i + 1] - y2a[i]) / 2.0;
    // 		double B = h * y2a[i] / 2.0 - (ya[i + 1] - ya[i]) / h;
    // 		double C = -h * y2a[i] / 6.0;

    // 		// Solve quadratic equation: A*t² + B*t + C = 0 where t = (x-xa[i])/h
    // 		double discriminant = B * B - 4 * A * C;

    // 		if (std::abs(A) < 1e-10) {
    // 			// Linear case
    // 			if (std::abs(B) > 1e-10) {
    // 				double t = -C / B;
    // 				double x = xa[i] + t * h;
    // 				if (x >= xMin && x <= xMax && x >= xa[i] && x <= xa[i + 1]) {
    // 					double y, dydx;
    // 					if (splint(xa, ya, y2a, n, x, &y, &dydx)) {
    // 						yMin = std::min(yMin, y);
    // 						yMax = std::max(yMax, y);
    // 					}
    // 				}
    // 			}
    // 		}
    // 		else if (discriminant >= 0) {
    // 			// Two possible roots
    // 			double t1 = (-B + sqrt(discriminant)) / (2 * A);
    // 			double t2 = (-B - sqrt(discriminant)) / (2 * A);

    // 			double x1 = xa[i] + t1 * h;
    // 			double x2 = xa[i] + t2 * h;

    // 			// Check if critical points are in this interval and the overall range
    // 			if (x1 >= xa[i] && x1 <= xa[i + 1] && x1 >= xMin && x1 <= xMax) {
    // 				double y, dydx;
    // 				if (splint(xa, ya, y2a, n, x1, &y, &dydx)) {
    // 					yMin = std::min(yMin, y);
    // 					yMax = std::max(yMax, y);
    // 				}
    // 			}

    // 			if (x2 >= xa[i] && x2 <= xa[i + 1] && x2 >= xMin && x2 <= xMax) {
    // 				double y, dydx;
    // 				if (splint(xa, ya, y2a, n, x2, &y, &dydx)) {
    // 					yMin = std::min(yMin, y);
    // 					yMax = std::max(yMax, y);
    // 				}
    // 			}
    // 		}
    // 	}

    // 	// Also check endpoints and knot points within range
    // 	for (int i = 0; i < n; i++) {
    // 		if (xa[i] >= xMin && xa[i] <= xMax) {
    // 			double y, dydx;
    // 			if (splint(xa, ya, y2a, n, xa[i], &y, &dydx)) {
    // 				yMin = std::min(yMin, y);
    // 				yMax = std::max(yMax, y);
    // 			}
    // 		}
    // 	}

    // 	// Check range endpoints if they're not knot points
    // 	double y, dydx;
    // 	if (splint(xa, ya, y2a, n, xMin, &y, &dydx)) {
    // 		yMin = std::min(yMin, y);
    // 		yMax = std::max(yMax, y);
    // 	}

    // 	if (splint(xa, ya, y2a, n, xMax, &y, &dydx)) {
    // 		yMin = std::min(yMin, y);
    // 		yMax = std::max(yMax, y);
    // 	}

    // 	// Check x = 0;
    // 	if (splint(xa, ya, y2a, n, 0, &y, &dydx)) {
    // 		yMin = std::min(yMin, y);
    // 		yMax = std::max(yMax, y);
    // 	}

    // 	return yMin != INFINITY && yMax != -INFINITY;
    // }

    // void process_cubic_spline_bounds(TElement* st_element, float x_minmax[2], float y_minmax[2],
    // 	float(&z_minmax)[2])
    // {
    // 	double z, x, y;

    // 	float& z_min = z_minmax[0];
    // 	float& z_max = z_minmax[1];

    // 	z_min = std::numeric_limits<double>::infinity();
    // 	z_max = -std::numeric_limits<double>::infinity();

    // 	float x_min_abs = get_absolute_minmax(x_minmax, 2, false);
    // 	float x_max_abs = get_absolute_minmax(x_minmax, 2, true);
    // 	float y_min_abs = get_absolute_minmax(y_minmax, 2, false);
    // 	float y_max_abs = get_absolute_minmax(y_minmax, 2, true);

    // 	double Rho_min = sqrt(x_min_abs * x_min_abs + y_min_abs * y_min_abs);
    // 	double Rho_max = sqrt(x_max_abs * x_max_abs + y_max_abs * y_max_abs);
    // 	double z_min_test;
    // 	double z_max_test;

    // 	find_spline_extrema(st_element->CubicSplineXData,
    // 		st_element->CubicSplineYData,
    // 		st_element->CubicSplineY2Data,
    // 		Rho_min, Rho_max, z_min_test, z_max_test);

    // 	/*for (int xi = 0; xi < 3; ++xi)
    // 	{
    // 		for (int yi = 0; yi < 3; ++yi)
    // 		{
    // 			x = x_range[xi];
    // 			y = y_range[yi];
    // 			double Rho = sqrt(x * x + y * y);
    // 			double dummy;
    // 			splint(st_element->CubicSplineXData,
    // 				st_element->CubicSplineYData,
    // 				st_element->CubicSplineY2Data,
    // 				st_element->CubicSplineXData.size(),
    // 				Rho, &z, &dummy);

    // 			if (z < z_min) z_min = z;
    // 			if (z > z_max) z_max = z;
    // 		}
    // 	}*/

    // 	z_min = z_min_test;
    // 	z_max = z_max_test;
    // }

    // void process_FE_bounds(TElement* st_element, float x_minmax[2], float y_minmax[2],
    // 	float(&z_minmax)[2])
    // {
    // 	float& z_min = z_minmax[0];
    // 	float& z_max = z_minmax[1];

    // 	z_min = std::numeric_limits<double>::infinity();
    // 	z_max = -std::numeric_limits<double>::infinity();

    // 	const MatDoub& xyz_nodes = st_element->FEData.nodes;

    // 	for (const std::vector<double>& fe_node : xyz_nodes)
    // 	{
    // 		if (fe_node[2] > z_max)
    // 			z_max = fe_node[2];
    // 		if (fe_node[2] < z_min)
    // 			z_min = fe_node[2];
    // 	}
    // }

    void transform_to_global(const float coord_element[3],
                             const tstage_ptr st_stage,
                             const TElement *st_element,
                             float (&coord_global)[3])
    {
        float PosDumStage[3];
        float coord_stage[3];
        MatrixVectorMult_generic(st_element->RLocToRef, coord_element, PosDumStage);
        for (int i = 0; i < 3; i++)
            coord_stage[i] = PosDumStage[i] + st_element->Origin[i];

        float PosDumGlob[3];
        MatrixVectorMult_generic(st_stage->RLocToRef, coord_stage, PosDumGlob);
        for (int i = 0; i < 3; i++)
            coord_global[i] = PosDumGlob[i] + st_stage->Origin[i];
        return;
    }

    void transform_bounds(const float min_coord_element[3],
                          const float max_coord_element[3],
                          const tstage_ptr st_stage,
                          const TElement *st_element,
                          float (&min_coord_global)[3],
                          float (&max_coord_global)[3])
    {
        // Transform min and max bounding box from element coordinates to global
        float corners_element[8][3] =
            {
                {min_coord_element[0], min_coord_element[1], min_coord_element[2]},
                {min_coord_element[0], min_coord_element[1], max_coord_element[2]},
                {min_coord_element[0], max_coord_element[1], min_coord_element[2]},
                {min_coord_element[0], max_coord_element[1], max_coord_element[2]},
                {max_coord_element[0], min_coord_element[1], min_coord_element[2]},
                {max_coord_element[0], min_coord_element[1], max_coord_element[2]},
                {max_coord_element[0], max_coord_element[1], min_coord_element[2]},
                {max_coord_element[0], max_coord_element[1], max_coord_element[2]}};

        // Convert corners to global coordinates
        float corners_global[8][3];
        for (int i = 0; i < 8; i++)
            transform_to_global(corners_element[i], st_stage, st_element, corners_global[i]);

        // Find min and max xyz
        min_coord_global[0] = corners_global[0][0];
        min_coord_global[1] = corners_global[0][1];
        min_coord_global[2] = corners_global[0][2];
        max_coord_global[0] = corners_global[0][0];
        max_coord_global[1] = corners_global[0][1];
        max_coord_global[2] = corners_global[0][2];
        for (int i = 1; i < 8; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                float val = corners_global[i][j];
                if (val < min_coord_global[j])
                    min_coord_global[j] = val;
                if (val > max_coord_global[j])
                    max_coord_global[j] = val;
            }
        }
    }

    bool get_bounds(const TElement *st_element,
                    float (&min_coord_global)[3],
                    float (&max_coord_global)[3])
    {
        // Get stage
        tstage_ptr st_stage = st_element->parent_stage;

        // Define element coord bounds
        float min_coord_element[3] = {0.f, 0.f, 0.f};
        float max_coord_element[3] = {0.f, 0.f, 0.f};

        double x_minmax[2] = {0.0, 0.0};
        double y_minmax[2] = {0.0, 0.0};
        double z_minmax[2] = {0.0, 0.0};

        st_element->aperture->bounding_box(x_minmax[0],
                                           x_minmax[1],
                                           y_minmax[0],
                                           y_minmax[1]);

        st_element->surface->bounding_box(x_minmax,
                                          y_minmax,
                                          z_minmax[0],
                                          z_minmax[1]);

        // Expand bounding boxes slightly to account for float precision
        const float expand = 1e-3f;
        x_minmax[0] -= expand;
        x_minmax[1] += expand;
        y_minmax[0] -= expand;
        y_minmax[1] += expand;
        z_minmax[0] -= expand;
        z_minmax[1] += expand;

        // Assign points to min/max coordinate element arrays
        min_coord_element[0] = x_minmax[0];
        min_coord_element[1] = y_minmax[0];
        min_coord_element[2] = z_minmax[0];
        max_coord_element[0] = x_minmax[1];
        max_coord_element[1] = y_minmax[1];
        max_coord_element[2] = z_minmax[1];

        // Convert local element bounds, to global xyz
        transform_bounds(min_coord_element, max_coord_element,
                         st_stage, st_element,
                         min_coord_global, max_coord_global);

        return true;
    }

} // namespace SolTrace::EmbreeRunner
