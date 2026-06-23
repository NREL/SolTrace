#pragma once

#include <QObject>
#include <QAbstractListModel>
#include "utilities/qt_helpers.h"

namespace SolTrace::GUI::App {
/**
 * @class StatusComponent
 * @brief Tracks the lifecycle state of an application module.
 *
 * Every app module owns a StatusComponent to communicate its readiness
 * to the QML layer. Status transitions follow a defined lifecycle:
 *
 * @code
 *   Unset → Loading → Incomplete → Ready → Complete → Stale → Error
 * @endcode
 *
 * - Unset:      Initial state, no data loaded
 * - Loading:    Actively reading from disk or backend
 * - Incomplete: Data loaded but required fields are missing
 * - Ready:      Data loaded and validated, ready for use
 * - Complete:   Fully configured and simulation-ready
 * - Stale:      Valid but potentially outdated (backend state changed)
 * - Error:      Validation failed or I/O error occurred
 *
 * QML binds to status to drive loading indicators, error banners,
 * and workflow gating (e.g. disabling "Run" if any module is Incomplete).
 */
class StatusComponent : public QObject {
    Q_OBJECT
public:
    explicit StatusComponent(QObject* parent = nullptr);

    enum class Status {
        Incomplete,
        Loading,
        Ready,
        Error,
        Stale,
        Complete,
        Unset
    };

    Q_ENUM(Status)

    Q_WRITABLE_PROPERTY(Status, status, Status::Unset)
    Q_WRITABLE_PROPERTY(QString, message, "")

    Q_INVOKABLE void mark_incomplete(const QString& message = "");
    Q_INVOKABLE void mark_loading(const QString& message = "");
    Q_INVOKABLE void mark_ready(const QString& message = "Ready");
    Q_INVOKABLE void mark_error(const QString& message);
    Q_INVOKABLE void mark_stale(const QString& message = "");
    Q_INVOKABLE void mark_complete(const QString& message = "");
    Q_INVOKABLE void mark_unset(const QString& message = "");

};

/**
 * @class PresetComponentBase
 * @brief QML-facing base class for preset management.
 *
 * Exposes preset operations to QML via Q_INVOKABLE methods.
 * Inherits QAbstractListModel so QML can display the list of available
 * presets directly (name + description per row).
 *
 * Actual file I/O is delegated to the backend's file utilities service.
 * This class is the thin QML adapter — it does not own persistence logic.
 *
 * Roles:
 * - NameRole:        Display name of the preset
 * - DescriptionRole: Optional description string
 */
class PresetComponentBase : public QAbstractListModel {
    Q_OBJECT
public:

    enum Roles { NameRole = Qt::UserRole + 1, DescriptionRole };
    QHash<int, QByteArray> roleNames() const override;

    Q_INVOKABLE bool create(const QString& name, const QString& description);
    Q_INVOKABLE bool load(const QString& name);
    Q_INVOKABLE bool modify(const QString& name);
    Q_INVOKABLE bool save(const QString& name, const QString& description = "");
    Q_INVOKABLE bool remove(const QString& name);
    Q_INVOKABLE bool import_preset(const QString& filepath);
    Q_INVOKABLE bool export_preset(const QString& name,
                                   const QString& filepath);

protected:
    virtual bool create_impl(const QString& name, const QString& description) = 0;
    virtual bool load_impl(const QString& name) = 0;
    virtual bool modify_impl(const QString& name) = 0;
    virtual bool save_impl(const QString& name, const QString& description = "") = 0;
    virtual bool remove_impl(const QString& name) = 0;
    virtual bool import_preset_impl(const QString& filepath) = 0;
    virtual bool export_preset_impl(const QString& name, const QString& filepath) = 0;
};

/**
 * @class PresetComponent
 * @brief Typed preset storage and active selection management.
 *
 * Template subclass of PresetComponentBase that adds typed storage
 * and tracks the currently active preset.
 *
 * @tparam T The domain type being managed (e.g. SunDefinition,
 * DirectionalSunPosition)
 *
 * active_preset reflects the currently loaded configuration. QML binds
 * to this to display and edit the active parameters. Changing active_preset
 * does not automatically save — the user must explicitly call save().
 *
 * m_keys preserves insertion order for stable row indexing in the list model.
 */

template <typename T>
struct Preset {
    QString description;
    QString file_path; // empty = never saved
    T* data;
    bool modified = false;
};

template <typename T>
class PresetComponent : public PresetComponentBase {

public:
    QPointer<T> active_preset();
    QPointer<T> set_active_preset(QString key);

    int rowCount(const QModelIndex &parent) const override
    {
        if (parent.isValid()) return 0;
        return m_presets.count();
    }

    QVariant data(const QModelIndex &index, int role) const override
    {
        if (!index.isValid()) return {};
        if (index.row() > m_keys.count()) return {};
        QString key = m_keys[index.row()];
        switch (role) {
        case Roles::NameRole: return key;
        case Roles::DescriptionRole: return m_presets[key].description;
        default: return {};
        }
    }

protected:
    bool create_impl(const QString& name, const QString& description) override {
        QString key_name = name;

        if (m_presets.contains(key_name)) {
            for (int i = 1; ; i++) {
                QString str = name + " (" + QString::number(i) + ")";
                if (!m_presets.contains(str)) {
                    key_name = str;
                    break;
                }
            }
        }

        m_presets.insert(key_name, Preset<T>{.description=description, .data = new T()});
        m_keys.append(key_name);
        m_active = key_name;

        return true;
    }

    bool modify_impl(const QString& name) override {
        if (!m_presets.contains(name)) return false;
        m_presets[name].modified = true;
        return true;
    }

    bool load_impl(const QString& name) override {
        if (!m_presets.contains(name)) return false;
        m_active = name;
        return true;
    }

    bool save_impl(const QString& name, const QString& description = "") override {
        if (!m_presets.contains(name)) return false;
        // stub
        return true;
    }

    bool remove_impl(const QString& name) override {
        return m_presets.remove(name);
    }

    bool import_preset_impl(const QString& filepath) override {
        // stub
        return false;
    }

    bool export_preset_impl(const QString& name, const QString& filepath) override {
        if (!m_presets.contains(name)) return false;
        // stub
        return false;
    }

private:
    QHash<QString,Preset<T>> m_presets;
    QList<QString>     m_keys;
    QString m_active;
};

} // namespace SolTrace::GUI::App
