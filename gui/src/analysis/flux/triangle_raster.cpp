#include "triangle_raster.h"

#include <algorithm>

namespace SolTrace::GUI::Analysis {

static double edge_function(QPointF const& a,
                            QPointF const& b,
                            QPointF const& p) {
    return (p.x() - a.x()) * (b.y() - a.y()) -
           (p.y() - a.y()) * (b.x() - a.x());
}

static bool is_top_left_edge(QPointF const& a, QPointF const& b) {
    auto const dx = b.x() - a.x();
    auto const dy = b.y() - a.y();
    return dy > 0.0 || (dy == 0.0 && dx < 0.0);
}

bool contains_with_edge_ownership(QPointF const& p,
                                  QPointF const& a,
                                  QPointF const& b,
                                  QPointF const& c) {
    auto aa = a;
    auto bb = b;
    auto cc = c;

    if (edge_function(aa, bb, cc) < 0.0) {
        std::swap(bb, cc);
    }

    auto accepts_edge = [](double edge_value,
                           QPointF const& edge_a,
                           QPointF const& edge_b) {
        constexpr double edge_epsilon = 1e-10;

        if (edge_value > edge_epsilon) return true;
        if (edge_value < -edge_epsilon) return false;

        return is_top_left_edge(edge_a, edge_b);
    };

    return accepts_edge(edge_function(aa, bb, p), aa, bb) &&
           accepts_edge(edge_function(bb, cc, p), bb, cc) &&
           accepts_edge(edge_function(cc, aa, p), cc, aa);
}

} // namespace SolTrace::GUI::Analysis
