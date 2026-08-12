#pragma once

#include "database/database.h"

#include <QMetaObject>
#include <QObject>
#include <QPointer>
#include <QVector>

namespace db {

/// Helper base for QObject adapters that observe a mutable Database.
///
/// Calling observe() disconnects previous database signal connections and lets
/// the derived class attach to the new database in
/// set_new_database_connections().
class DatabaseObserver {
    QPointer<Database>               m_database;
    QVector<QMetaObject::Connection> m_database_conns;

protected:
    void observe(Database* ptr) {
        if (ptr == m_database) return;
        if (m_database) {
            for (auto const& c : std::as_const(m_database_conns)) {
                QObject::disconnect(c);
            }
            m_database_conns.clear();
        }
        m_database = ptr;
        if (ptr) set_new_database_connections(ptr);
    }

    void add_connection(QMetaObject::Connection c) {
        m_database_conns.push_back(c);
    }

    Database*       database() { return m_database; }
    Database const* database() const { return m_database; }

    virtual void set_new_database_connections(Database* ptr) = 0;

    template <class F>
    void with_db(F&& f) {
        if (m_database) { f(m_database); }
    }

public:
    DatabaseObserver()          = default;
    virtual ~DatabaseObserver() = default;

    DatabaseObserver(DatabaseObserver const&)            = delete;
    DatabaseObserver& operator=(DatabaseObserver const&) = delete;
    DatabaseObserver(DatabaseObserver&&)                 = delete;
    DatabaseObserver& operator=(DatabaseObserver&&)      = delete;
};

/// Helper base for QObject adapters that observe a const Database.
class ConstDatabaseObserver {
    QPointer<Database const>         m_database;
    QVector<QMetaObject::Connection> m_database_conns;

protected:
    void observe(Database const* ptr) {
        if (ptr == m_database) return;
        if (m_database) {
            for (auto const& c : std::as_const(m_database_conns)) {
                QObject::disconnect(c);
            }
            m_database_conns.clear();
        }
        m_database = ptr;
        if (ptr) set_new_database_connections(ptr);
    }

    void add_connection(QMetaObject::Connection c) {
        m_database_conns.push_back(c);
    }

    Database const* database() const { return m_database; }

    virtual void set_new_database_connections(Database const* ptr) = 0;

    template <class F>
    void with_db(F&& f) {
        if (m_database) { f(m_database); }
    }

public:
    ConstDatabaseObserver()          = default;
    virtual ~ConstDatabaseObserver() = default;

    ConstDatabaseObserver(ConstDatabaseObserver const&)            = delete;
    ConstDatabaseObserver& operator=(ConstDatabaseObserver const&) = delete;
    ConstDatabaseObserver(ConstDatabaseObserver&&)                 = delete;
    ConstDatabaseObserver& operator=(ConstDatabaseObserver&&)      = delete;
};

} // namespace db
