#pragma once

#include "script/script_db_interface.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

#include <QObject>
#include <QQmlEngine>
#include <QStringList>

namespace SolTrace::GUI::Script {

/// One user-editable argument declared in a script header.
///
/// Script headers are parsed into these records and exposed to QML through
/// ScriptPropertyModel. The `value` field is editable; the remaining fields
/// describe how the UI should present and validate that value.
struct ScriptProperty {
    QString name;
    QString identifier;
    QString type;
    QString extra;
    bool    min_bounded   = false;
    bool    max_bounded   = false;
    bool    max_inclusive = false;
    double  min           = 0.0;
    double  max           = 0.0;
    QString unit;
    QString error;

    QString value;

    RECORD_META(ScriptProperty,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RO(identifier),
                SM_EXPOSE_RO(type),
                SM_EXPOSE_RO(extra),
                SM_EXPOSE_RO(min_bounded),
                SM_EXPOSE_RO(max_bounded),
                SM_EXPOSE_RO(max_inclusive),
                SM_EXPOSE_RO(min),
                SM_EXPOSE_RO(max),
                SM_EXPOSE_RO(unit),
                SM_EXPOSE_RO(error),
                SM_EXPOSE_RW(value));
};

class ScriptPropertyModel : public StructTableModel<ScriptProperty> {
    Q_OBJECT

public:
    explicit ScriptPropertyModel(QObject* parent = nullptr);
};

/// User-authored script plus parsed metadata and execution state.
///
/// Script source must start with a comment header block. No whitespace or code
/// may appear before it. Both block comments and consecutive line comments are
/// accepted:
///
///     /*
///     TITLE Example Script
///     DESC Creates a few entities.
///     DESC Additional DESC lines are concatenated with newlines.
///     PROPERTY mirror_count integer 12 1..=100
///     PROPERTY radius real 10.0 0.1..
///     PROPERTY label string Demo
///     PROPERTY direction vec3 [0,0,1] unit
///     *\/
///     const e = db.create()
///     return e
///
///     // TITLE Example Script
///     // DESC Equivalent line-comment form.
///     // PROPERTY count integer 4 1..=10
///     return db.create()
///
/// Header directives:
/// - TITLE text
///   Required. Must be the first non-empty directive in the header. The
///   remaining text becomes Script::title.
/// - DESC text
///   Optional and repeatable. Each line is appended to Script::description.
/// - PROPERTY name type initial [extra...]
///   Optional and repeatable. Properties are parsed in header order, exposed in
///   Script::properties, and passed as positional arguments to the evaluated
///   JavaScript function during run(). The initial value is stored as the
///   editable ScriptProperty::value shown in the UI.
///
/// Property types:
/// - integer initial range
///   Integer argument. The range follows the initial value and uses Rust-like
///   syntax: `..0`, `1..`, `1..10`, or `1..=100`. A range is currently
///   required.
/// - real initial range
///   Floating-point argument with the same range syntax as integer.
/// - string
///   Unrestricted string argument. The initial value is one token; no extra
///   arguments are currently used.
/// - vec3 initial [unit]
///   Three-component vector argument. The optional `unit` extra follows the
///   initial value and marks that the vector should be treated as a direction.
///   The current runner accepts values like `1 0 0`, `[1, 0, 0]`, or a single
///   scalar such as `1`, which expands to `[1, 1, 1]`. Header initial values
///   should be written without spaces, for example `[0,0,1]`.
///
/// After the header, the script body is wrapped in a generated JavaScript
/// function. Script::run() converts property values according to the
/// declarations above and calls the generated function with those arguments.
/// PROPERTY identifiers are available as normal JavaScript parameters, and the
/// database script API is available as `db`.
class Script : public QObject {
    Q_OBJECT

    QPointer<ScriptDBInterface> m_interface;
    QPointer<db::Database>      m_database;

    Q_WRITABLE_PROPERTY(QString, code, {});
    Q_READONLY_PROPERTY(QString, title);
    Q_READONLY_PROPERTY(QString, description);
    Q_READONLY_PROPERTY(bool, valid);
    Q_READONLY_PROPERTY(QStringList, parse_errors);
    Q_READONLY_PROPERTY(QStringList, run_errors);
    QOBJECT_READONLY_PROPERTY(ScriptPropertyModel, properties);

    Q_READONLY_PROPERTY(QStringList, builtin_scripts);

public:
    explicit Script(QObject* parent = nullptr);

    void set_database(db::Database*);

public slots:
    bool parse();

    void run();
    void notify_error(QString message);

signals:
    void notify(ANotification);
};

} // namespace SolTrace::GUI::Script
