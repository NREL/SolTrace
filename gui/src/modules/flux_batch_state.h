#pragma once

#include "data/entity.h"

#include <QSet>
#include <QtTypes>

namespace SolTrace::GUI::App {

class FluxBatchState {
    qint64           m_max = 0;
    QSet<SolTrace::GUI::Data::Entity> m_pending;

public:
    void reset(qint64 max = 0);
    void add_pending(SolTrace::GUI::Data::Entity entity);
    void start_from_pending();
    bool is_pending() const;
    bool empty() const;
    qint64 max() const;
    qint64 remaining() const;
    qint64 completed() const;
    bool mark_done(SolTrace::GUI::Data::Entity entity);
};

} // namespace SolTrace::GUI::App
