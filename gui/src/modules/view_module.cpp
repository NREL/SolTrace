#include "view_module.h"

namespace SolTrace::GUI::App {

ViewModule::ViewModule(QObject* parent)
    : QObject { parent },
      m_left_panel(new SplitPanelData(this)),
      m_right_panel(new SplitPanelData(this)),
      m_full_panel(new FullPanelData(this)),
      m_sim(new SimulationViewState(this)) { }

void SplitPanelData::update_size() {
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

void SplitPanelData::save_visibility() {
    set_saved_visible(m_visible);
}

void SplitPanelData::restore_visibility() {
    set_visible(m_saved_visible);
}

void SplitPanelData::show() {
    set_visible(true);
}

void SplitPanelData::hide() {
    set_visible(false);
}

SplitPanelData::SplitPanelData(QObject* parent) : QObject(parent) {
    connect(this,
            &SplitPanelData::width_changed,
            this,
            &SplitPanelData::update_size);
    update_size();
}

QVector<int> SplitPanelData::sizes() const {
    return m_sizes;
}

QVector<int> SplitPanelData::thresholds() const {
    return m_thresholds;
}

bool SplitPanelData::is_small() {
    return m_size == Small;
}
bool SplitPanelData::is_normal() {
    return m_size == Normal;
}
bool SplitPanelData::is_wide() {
    return m_size == Wide;
}

bool ViewModule::shrink_panel(const QVector<int>&       sizes,
                              QPointer<SplitPanelData>& p) {
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
    QPointer<SplitPanelData> expanding =
        expanding_right_panel ? m_right_panel : m_left_panel;
    QPointer<SplitPanelData> collapsing =
        expanding_right_panel ? m_left_panel : m_right_panel;

    const auto& sizes   = expanding->sizes();
    const int   small_w = sizes[SplitPanelData::Small];
    const int   full_w  = sizes[SplitPanelData::Full];

    // Short-circuit cases:
    // 1) Full mode: expands self and hides other
    if (expanding->size() == SplitPanelData::Full && !resizing_window) {
        expanding->set_width(full_w);
        expanding->show();
        collapsing->hide();
        return;
    }

    // 2) Other panel is full — knock it down so we can fit
    if (collapsing->visible() && collapsing->size() == SplitPanelData::Full &&
        !resizing_window) {
        collapsing->set_width(sizes[SplitPanelData::Wide]);
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
            collapsing->hide();
        }
        // Last resort: shrink the expanding panel itself
        else if (expanding->width() > small_w && resizing_window) {
            if (!shrink_panel(sizes, expanding)) break;
        } else {
            break;
        }
    }
}

void ViewModule::open_full_panel() {
    m_left_panel->save_visibility();
    m_right_panel->save_visibility();
    m_left_panel->hide();
    m_right_panel->hide();
    m_full_panel->show();
}

void ViewModule::close_full_panel(int available_width) {
    m_left_panel->restore_visibility();
    m_right_panel->restore_visibility();
    fit_panels(available_width);
    m_full_panel->hide();
}

void ViewModule::toggle_full_panel(int available_width) {
    if (m_full_panel->visible()) {
        close_full_panel(available_width);
        return;
    } else {
        open_full_panel();
    }
}

void ViewModule::open_left_panel(int available_width) {
    if (m_full_panel->visible()) {
        close_full_panel(available_width);
        return;
    }
    m_left_panel->show();
    fit_panels(available_width);
}

void ViewModule::close_left_panel() {
    m_left_panel->hide();
}

void ViewModule::toggle_left_panel(int available_width) {
    if (m_left_panel->visible()) {
        close_left_panel();
        return;
    }
    open_left_panel(available_width);
}

void ViewModule::open_right_panel(int available_width) {
    if (m_full_panel->visible()) {
        close_full_panel(available_width);
        return;
    }
    m_right_panel->show();
    fit_panels(available_width, true);
}

void ViewModule::close_right_panel() {
    m_right_panel->hide();
}

void ViewModule::toggle_right_panel(int available_width) {
    if (m_right_panel->visible()) {
        close_right_panel();
        return;
    }
    open_right_panel(available_width);
}


FullPanelData::FullPanelData(QObject* parent) : QObject(parent) { }

void FullPanelData::show() {
    set_visible(true);
}

void FullPanelData::hide() {
    set_visible(false);
}


} // namespace SolTrace::GUI::App
