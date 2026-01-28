
#include "common.hpp"

#include <cmath>

#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

bool is_identical(const Vector3d &x, const Vector3d &y)
{
    return (
        x.data[0] == y.data[0] &&
        x.data[1] == y.data[1] &&
        x.data[2] == y.data[2]);
}

bool is_identical(const Vector3d &x, const Vector3d &y, double tol)
{
    return (
        fabs(x.data[0] - y.data[0]) <= tol &&
        fabs(x.data[1] - y.data[1]) <= tol &&
        fabs(x.data[2] - y.data[2]) <= tol);
}

bool is_identical(const Matrix3d &A, const Matrix3d &B)
{
    bool all_identical = true;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            all_identical &= A.data[i][j] == B.data[i][j];
        }
    }
    return all_identical;
}

element_ptr make_configured_element()
{
    element_ptr el = SolTrace::Data::make_element<SolTrace::Data::SingleElement>();
    el->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Circle>(2.0));
    el->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    return el;
}

int_fast64_t count_element_event(const SimulationResult &res, element_id el, RayEvent rev)
{
    int_fast64_t count = 0;

    for (auto ray_idx = 0;
         ray_idx < res.get_number_of_records();
         ++ray_idx)
    {
        auto rr = res[ray_idx];
        for (auto event_idx = 0;
             event_idx < rr->get_number_of_interactions();
             ++event_idx)
        {
            if (rr->get_event(event_idx) == rev &&
                rr->get_element(event_idx) == el)
            {
                ++count;
            }
        }
    }

    return count;
}
