#include "view_module.h"

namespace SolTrace::GUI::App {

ViewModule::ViewModule(QObject* parent)
    : QObject { parent },
      m_left_panel(new PanelData()),
      m_right_panel(new PanelData()),
      m_settings_panel(new PanelData()),
      m_sim(new SimulationViewState()) { }

void PanelData::update_size() {
    PanelSize new_size;
    if (m_width < m_thresholds[0]) new_size = Small;
    else if (m_width < m_thresholds[1])
        new_size = Normal;
    else if (m_width < m_thresholds[2])
        new_size = Wide;
    else
        new_size = Full;

    set_size(new_size);
}

PanelData::PanelData(QObject* parent) : QObject(parent) {
    connect(this, &PanelData::width_changed, this, &PanelData::update_size);
    update_size();
}

QVector<int> PanelData::sizes() const {
    return m_sizes;
}

QVector<int> PanelData::thresholds() const {
    return m_thresholds;
}

bool PanelData::is_small() {
    return m_size == Small;
}
bool PanelData::is_normal() {
    return m_size == Normal;
}
bool PanelData::is_wide() {
    return m_size == Wide;
}

bool ViewModule::shrink_panel(const QVector<int>&  sizes,
                              QPointer<PanelData>& p) {
    const int w = p->width();
    // Find the largest preset that is strictly less than current width.
    for (int i = sizes.size() - 1; i >= 0; --i) {
        if (sizes[i] < w) {
            p->set_width(sizes[i]);
            return true;
        }
    }
    return false;
}


void ViewModule::fit_panels(int  available_width,
                            bool expanding_right_panel,
                            bool resizing_window,
                            int  margin) {
    QPointer<PanelData> expanding =
        expanding_right_panel ? m_right_panel : m_left_panel;
    QPointer<PanelData> collapsing =
        expanding_right_panel ? m_left_panel : m_right_panel;

    const auto& sizes   = expanding->sizes();
    const int   small_w = sizes[PanelData::Small];
    const int   full_w  = sizes[PanelData::Full];

    // Short-circuit cases:
    // 1) Full mode: expands self and hides other
    if (expanding->size() == PanelData::Full && !resizing_window) {
        expanding->set_width(full_w);
        expanding->set_visible(true);
        collapsing->set_visible(false);
        return;
    }

    // 2) Other panel is full — knock it down so we can fit
    if (collapsing->visible() && collapsing->size() == PanelData::Full &&
        !resizing_window) {
        collapsing->set_width(sizes[PanelData::Wide]);
    }

    while (true) {
        int left_w  = m_left_panel->visible() ? m_left_panel->width() : 0;
        int right_w = m_right_panel->visible() ? m_right_panel->width() : 0;

        if (left_w + right_w + margin <= available_width) break;

        // First try shrinking the non-expanding panel
        if (collapsing->visible() && collapsing->width() > small_w) {
            if (!shrink_panel(sizes, collapsing)) break;
        }
        // If it can't shrink further, hide it
        else if (collapsing->visible() && !resizing_window) {
            collapsing->set_visible(false);
        }
        // Last resort: shrink the expanding panel itself
        else if (expanding->width() > small_w && resizing_window) {
            if (!shrink_panel(sizes, expanding)) break;
        } else {
            break;
        }
    }
}

} // namespace SolTrace::GUI::App
