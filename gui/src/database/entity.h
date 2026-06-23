#pragma once

#include <entt/entt.hpp>

#include <QDebug>
#include <QString>
#include <qqmlintegration.h>

namespace db {


/// Wrapper type for Qt/QML interface
struct Entity {
    Q_GADGET
    QML_VALUE_TYPE(db_entity);
    // Q_PROPERTY(qint64 value MEMBER value);

public:
    Entity() = default;
    Entity(entt::entity e) : value(e) { }

    entt::entity value = entt::null;

    std::strong_ordering operator<=>(Entity const& other) const = default;

    operator entt::entity() const { return value; }

    Q_INVOKABLE bool is_valid() const { return value != entt::null; }

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
struct hash<db::Entity> {
    std::size_t operator()(db::Entity e) const noexcept {
        return std::hash<entt::entity>()(e.value);
    }
};
} // namespace std
