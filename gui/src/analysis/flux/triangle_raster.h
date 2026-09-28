#pragma once

#include <QPointF>

namespace SolTrace::GUI::Analysis {

bool contains_with_edge_ownership(QPointF const& p,
                                  QPointF const& a,
                                  QPointF const& b,
                                  QPointF const& c);

} // namespace SolTrace::GUI::Analysis
