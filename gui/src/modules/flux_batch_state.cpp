#include "flux_batch_state.h"

namespace SolTrace::GUI::App {

void FluxBatchState::reset(qint64 max) {
    m_max = max;
    m_pending.clear();
}

void FluxBatchState::add_pending(SolTrace::GUI::Data::Entity entity) {
    m_pending.insert(entity);
}

void FluxBatchState::start_from_pending() {
    m_max = m_pending.size();
}

bool FluxBatchState::is_pending() const {
    return !m_pending.empty();
}

bool FluxBatchState::empty() const {
    return m_pending.empty();
}

qint64 FluxBatchState::max() const {
    return m_max;
}

qint64 FluxBatchState::remaining() const {
    return m_pending.size();
}

qint64 FluxBatchState::completed() const {
    if (m_max <= 0) return 0;
    return m_max - remaining();
}

bool FluxBatchState::mark_done(SolTrace::GUI::Data::Entity entity) {
    if (m_max <= 0) return false;
    if (!m_pending.contains(entity)) return false;

    m_pending.remove(entity);
    return true;
}

} // namespace SolTrace::GUI::App
