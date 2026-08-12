#pragma once

#include <entt/entt.hpp>

#include <QDebug>
#include <QString>
#include <qqmlintegration.h>

namespace db {

/// Wrapper type for passing entt entities through Qt/QML.
///
/// Entity values are only meaningful inside the database that created them.
struct Entity {
    Q_GADGET
    QML_VALUE_TYPE(db_entity);

public:
    Entity() = default;

    /// Wrap an entt entity handle.
    Entity(entt::entity e) : value(e) { }

    entt::entity value = entt::null;

    std::strong_ordering operator<=>(Entity const& other) const = default;

    operator entt::entity() const { return value; }

    /// Whether this wrapper contains a non-null entity handle.
    Q_INVOKABLE bool is_valid() const { return value != entt::null; }

    /// Debug-friendly representation for QML and logs.
    Q_INVOKABLE QString debug_string() const {
        if (!is_valid()) return QStringLiteral("entity(null)");
        return QStringLiteral("entity(%1)").arg(entt::to_integral(value));
    }
};

inline QDebug operator<<(QDebug debug, Entity const& c) {
    QDebugStateSaver saver(debug);
    debug.nospace() << "(Entity " << entt::to_integral(c.value) << ")";

    return debug;
}

} // namespace db


namespace std {
template <>
/// Allow db::Entity keys in unordered containers.
struct hash<db::Entity> {
    std::size_t operator()(db::Entity e) const noexcept {
        return std::hash<entt::entity>()(e.value);
    }
};
} // namespace std
