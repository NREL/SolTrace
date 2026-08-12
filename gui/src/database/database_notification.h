#pragma once

#include "components.h"

#include <type_traits>
#include <utility>

#include <QDebug>
#include <QObject>
#include <QTimer>

#include <entt/entity/fwd.hpp>

namespace db {


/// Helper function: patch a component, skipping it if it does not exist.
/// Returns true if the patch occurred.
template <class Component, class Function>
bool try_patch(entt::registry& reg, entt::entity entity, Function&& f) {
    if (!reg.valid(entity)) return false;

    if (reg.all_of<Component>(entity)) {
        reg.patch<Component>(entity, f);
        return true;
    }

    return false;
}

/// Helper function: patch a component, creating it if it does not exist.
template <class Component, class Function>
void emplace_patch(entt::registry& reg, entt::entity entity, Function&& f) {
    if (!reg.all_of<Component>(entity)) {
        if constexpr (std::is_empty_v<Component>) {
            reg.emplace<Component>(entity);
        } else {
            reg.emplace<Component>(entity, Component { });
        }
    }

    reg.patch<Component>(entity, f);
}

/// Base class for all component notification, needed as Qt does not support
/// signals+slots on template classes.
class ComponentAPIBase : public QObject {
    Q_OBJECT
protected:
    entt::registry& m_host;

public:
    ComponentAPIBase(entt::registry& p) : m_host(p) { }
    virtual ~ComponentAPIBase() = default;

signals:
    /// An entity with this component has had that component changed
    void changed(entt::entity);
    /// An entity with this component has had that component removed
    void removed(entt::entity);
};

/// Notification endpoint for a given component.
///
/// Bridges EnTT construct/update/destroy callbacks into Qt signals.
template <class Component>
class ComponentAPI : public ComponentAPIBase {
    void change_callback(entt::registry& reg, entt::entity entity) {
        emit changed(entity);
    }

    void remove_callback(entt::registry& reg, entt::entity entity) {
        // We HAVE to do this, because entt's callback is "about to be removed"
        // not "removed"
        QTimer::singleShot(
            0, this, [this, entity]() { emit this->removed(entity); });
    }

public:
    ComponentAPI(entt::registry& p) : ComponentAPIBase(p) {
        p.on_construct<Component>()
            .template connect<&ComponentAPI::change_callback>(this);
        p.on_update<Component>()
            .template connect<&ComponentAPI::change_callback>(this);
        p.on_destroy<Component>()
            .template connect<&ComponentAPI::remove_callback>(this);
    }
    ~ComponentAPI() {
        m_host.template on_construct<Component>()
            .template disconnect<&ComponentAPI::change_callback>(this);
        m_host.template on_update<Component>()
            .template disconnect<&ComponentAPI::change_callback>(this);
        m_host.template on_destroy<Component>()
            .template disconnect<&ComponentAPI::remove_callback>(this);
    }

    /// A self pointer, useful for connect() calls.
    auto* self() { return this; }

    /// Get the content of this component on an entity. Returns null if the
    /// content is not available
    Component const* get(entt::entity entity) const {
        if (m_host.valid(entity)) {
            return m_host.template try_get<Component>(entity);
        }
        return nullptr;
    }

    /// Patch a component if it already exists on entity.
    template <class F>
    bool try_patch(entt::entity entity, F&& f) {
        if (!m_host.valid(entity)) return false;

        if (m_host.all_of<Component>(entity)) {
            m_host.patch<Component>(entity, f);
            return true;
        }

        return false;
    }

    /// Patch a component, creating it first if needed.
    template <class F>
    void emplace_patch(entt::entity entity, F&& f) {
        if (!m_host.all_of<Component>(entity)) {
            if constexpr (std::is_empty_v<Component>) {
                m_host.emplace<Component>(entity);
            } else {
                m_host.emplace<Component>(entity, Component { });
            }
        }

        m_host.patch<Component>(entity, f);
    }

    /// A const view of all entities with this component
    auto view() const { return m_host.view<Component const>(); }
};

/// Specialized component notification endpoint that supports editing.
template <class Component>
class ComponentAPIUpdate : public ComponentAPI<Component> {
public:
    ComponentAPIUpdate(entt::registry& p) : ComponentAPI<Component>(p) { }
    ~ComponentAPIUpdate() override = default;

public:
    /// Replace the component value on entity.
    void set(entt::entity entity, Component const& c) {
        if constexpr (std::is_empty_v<Component>) {
            this->m_host.template emplace_or_replace<Component>(entity);
        } else {
            this->m_host.template emplace_or_replace<Component>(entity, c);
        }
    }

    /// Patch the component value on entity, creating it first if needed.
    template <class Function>
    void patch(entt::entity entity, Function&& f) {
        db::emplace_patch<Component>(this->m_host, entity, f);
    }

    /// Remove the component from entity.
    void remove(entt::entity entity) {
        this->m_host.template remove<Component>(entity);
    }
};

/// Notification endpoint for a component that we treat as a unique resource.
template <class Component>
class SingletonComponentAPI : public ComponentAPIBase {
    entt::entity m_entity      = entt::null;
    bool         m_owns_entity = false;

    void verify_singleton() const {
        auto view = this->m_host.template view<Component const>();
        auto it   = view.begin();

        if (it == view.end()) return;

        const auto first = *it;
        ++it;

        if (it != view.end()) {
            qFatal("Singleton component exists on more than one entity");
        }

        if (m_entity != entt::null && m_entity != first) {
            qFatal("Singleton component entity tracking is inconsistent");
        }
    }

    void change_callback(entt::registry& reg, entt::entity entity) {
        if (m_entity != entt::null && m_entity != entity) {
            qFatal("Singleton component exists on more than one entity");
        }

        m_entity = entity;
        emit changed(entity);
    }

    void remove_callback(entt::registry& reg, entt::entity entity) {
        if (m_entity == entity) {
            m_entity      = entt::null;
            m_owns_entity = false;
        }

        QTimer::singleShot(
            0, this, [this, entity]() { emit this->removed(entity); });
    }

public:
    explicit SingletonComponentAPI(entt::registry& p) : ComponentAPIBase(p) {
        auto view = p.template view<Component const>();
        for (auto entity : view) {
            if (m_entity != entt::null) {
                qFatal("Singleton component exists on more than one entity");
            }
            m_entity      = entity;
            m_owns_entity = false;
        }

        p.on_construct<Component>()
            .template connect<&SingletonComponentAPI::change_callback>(this);
        p.on_update<Component>()
            .template connect<&SingletonComponentAPI::change_callback>(this);
        p.on_destroy<Component>()
            .template connect<&SingletonComponentAPI::remove_callback>(this);

        verify_singleton();
    }

    ~SingletonComponentAPI() override {
        this->m_host.template on_construct<Component>()
            .template disconnect<&SingletonComponentAPI::change_callback>(this);
        this->m_host.template on_update<Component>()
            .template disconnect<&SingletonComponentAPI::change_callback>(this);
        this->m_host.template on_destroy<Component>()
            .template disconnect<&SingletonComponentAPI::remove_callback>(this);
    }

    /// A self pointer, useful for connect() calls.
    auto* self() { return this; }

    /// Get the entity that owns this singleton
    entt::entity entity() const { return m_entity; }

    /// Does this singleton exist?
    bool exists() const { return m_entity != entt::null; }

    /// Get the component, or null if it does not exist
    Component const* get() const {
        if (m_entity == entt::null) return nullptr;
        if (!this->m_host.valid(m_entity)) return nullptr;
        return this->m_host.template try_get<Component>(m_entity);
    }

    /// Get the component, or null if it does not exist
    Component* get() {
        if (m_entity == entt::null) return nullptr;
        if (!this->m_host.valid(m_entity)) return nullptr;
        return this->m_host.template try_get<Component>(m_entity);
    }

    /// Get the component, panic if it does not exist
    Component const& require() const {
        auto ptr = get();
        if (!ptr) qFatal("Required singleton component is missing");
        return *ptr;
    }

    /// Get the component, panic if it does not exist
    Component& require() {
        auto ptr = get();
        if (!ptr) qFatal("Required singleton component is missing");
        return *ptr;
    }

    /// Set the component content, replacing it if it exists
    void set(Component const& component) {
        if (m_entity == entt::null) {
            m_entity      = this->m_host.create();
            m_owns_entity = true;
        }

        if constexpr (std::is_empty_v<Component>) {
            this->m_host.template emplace_or_replace<Component>(m_entity);
        } else {
            this->m_host.template emplace_or_replace<Component>(m_entity,
                                                                component);
        }
    }

    /// Run a patch function on the singleton, emplacing it with a default if it
    /// does not exist
    template <class Function>
    void patch(Function&& f) {
        if (m_entity == entt::null) {
            m_entity      = this->m_host.create();
            m_owns_entity = true;
        }
        db::emplace_patch<Component>(
            this->m_host, m_entity, std::forward<Function>(f));
    }

    /// Remove the singleton from the world.
    void remove() {
        if (m_entity == entt::null) return;

        const auto entity      = m_entity;
        const bool owns_entity = m_owns_entity;
        if (this->m_host.valid(entity)) {
            this->m_host.template remove<Component>(entity);
            if (owns_entity) { this->m_host.destroy(entity); }
        }
        m_entity      = entt::null;
        m_owns_entity = false;
    }
};

} // namespace db
