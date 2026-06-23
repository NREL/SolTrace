#pragma once


#include "opticaleditor.h"
#include "utilities/qt_helpers.h"

#include <QObject>

namespace db {

class MaterialEditor : public QObject, public DatabaseObserver {
    Q_OBJECT

    entt::entity m_current_group = entt::null;

    void set_new_database_connections(Database* ptr) override;


    QOBJECT_WRITABLE_PROPERTY(OpticalPropertiesObject, back_editor);
    QOBJECT_WRITABLE_PROPERTY(OpticalPropertiesObject, front_editor);

    // UX Helpers
    QOBJECT_READONLY_PROPERTY(QStringListModel, interaction_type_model);
    QOBJECT_READONLY_PROPERTY(QStringListModel, distribution_type_model);

private slots:
    void parameters_changed(entt::entity);

public:
    explicit MaterialEditor(QObject* parent = nullptr);
    ~MaterialEditor() override;

    void set(Database*, entt::entity group);

signals:

    void updated();
};

} // namespace db
