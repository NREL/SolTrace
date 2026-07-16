#include "script.h"
#include "script/schema_builder.h"

#include <QDebug>
#include <QDir>
#include <QDirIterator>
#include <QRegularExpression>

namespace SolTrace::GUI::Script {

namespace {

struct HeaderBlock {
    QString text;
    qsizetype body_start = 0;
    QString error;
};

QString normalize_header_line(QString line) {
    line = line.trimmed();
    if (line.startsWith('*')) { line = line.mid(1).trimmed(); }
    return line;
}

HeaderBlock first_comment_block(QString const& code) {
    if (code.startsWith("/*")) {
        auto end = code.indexOf("*/", 2);
        if (end < 0) {
            return { {}, 0, "Unterminated block comment header" };
        }
        return { code.mid(2, end - 2), end + 2, {} };
    }

    if (!code.startsWith("//")) {
        return { {}, 0, "Script must start with a comment header" };
    }

    QStringList lines;
    qsizetype   position = 0;

    while (position < code.size()) {
        auto line_end = code.indexOf('\n', position);
        if (line_end < 0) { line_end = code.size(); }

        auto line = code.mid(position, line_end - position);
        if (line.endsWith('\r')) { line.chop(1); }
        if (!line.startsWith("//")) { break; }
        lines << line.mid(2);
        position = line_end < code.size() ? line_end + 1 : line_end;
    }

    return { lines.join('\n'), position, {} };
}

bool parse_number(QString const& text, double& value) {
    bool ok = false;
    value   = text.toDouble(&ok);
    return ok;
}

void parse_range(QString const& range, ScriptProperty& property) {

    auto parts = range.split("..", Qt::KeepEmptyParts);

    if (parts.empty()) {
        property.error = "Invalid range";
        return;
    }

    auto lower = parts.value(0);
    auto upper = parts.value(1);

    if (upper.contains("=")) {
        upper                  = upper.replace("=", "");
        property.max_inclusive = true;
    }

    if (!lower.isEmpty()) {
        property.min_bounded = true;
        if (!parse_number(lower, property.min)) {
            property.error = "Invalid lower bound";
            return;
        }
    }

    if (!upper.isEmpty()) {
        property.max_bounded = true;
        if (!parse_number(upper, property.max)) {
            property.error = "Invalid upper bound";
            return;
        }
    }

    if (!property.min_bounded && !property.max_bounded) {
        property.error = "Range must include at least one bound";
        return;
    }

    if (property.min_bounded && property.max_bounded &&
        property.min > property.max) {
        property.error = "Lower bound exceeds upper bound";
    }
}

ScriptProperty parse_property(QString const& line) {
    auto parts = line.simplified().split(' ', Qt::SkipEmptyParts);

    ScriptProperty property;
    if (parts.size() < 4) {
        property.error = "PROPERTY requires name, type, and initial value";
        return property;
    }

    static QRegularExpression const identifier_pattern {
        QStringLiteral("^[A-Za-z_$][A-Za-z0-9_$]*$")
    };

    property.identifier = parts[1];
    if (!identifier_pattern.match(property.identifier).hasMatch()) {
        property.error =
            "PROPERTY name must be a valid JavaScript identifier";
        return property;
    }

    auto name = property.identifier.toLower();
    name.replace("_", " ");
    name[0] = name[0].toUpper();

    property.name = name;
    property.type = parts[2].toLower();
    property.value = parts[3];
    if (parts.size() > 4) { property.extra = parts.mid(4).join(' '); }

    if (property.name.isEmpty()) {
        property.error = "PROPERTY name is empty";
        return property;
    }

    if (property.type == "integer" || property.type == "real") {
        if (property.extra.isEmpty()) {
            property.error = "Numeric PROPERTY requires a range";
            return property;
        }
        parse_range(property.extra, property);
        return property;
    }

    if (property.type == "vec3") {
        if (!property.extra.isEmpty() && property.extra != "unit") {
            property.error = "vec3 PROPERTY only supports unit";
            return property;
        }
        property.unit = property.extra;
        return property;
    }

    if (property.type == "string") { return property; }

    property.error = "Unsupported PROPERTY type";
    return property;
}

QString script_function(QString const& code,
                        qsizetype      body_start,
                        QVector<ScriptProperty> const& properties) {
    QStringList arguments;
    for (auto const& property : properties) {
        arguments << property.identifier;
    }

    qsizetype body_start_line = 1;
    for (qsizetype i = 0; i < body_start && i < code.size(); ++i) {
        if (code.at(i) == QLatin1Char('\n')) { ++body_start_line; }
    }

    auto body = code.mid(body_start);
    return QStringLiteral("(function(%1) {%2%3\n})")
        .arg(arguments.join(QStringLiteral(", ")),
             QString(body_start_line - 1, QLatin1Char('\n')),
             body);
}

bool looks_like_legacy_function_script(QString const& code,
                                       qsizetype      body_start) {
    auto body = code.mid(body_start).trimmed();
    return body.startsWith(QStringLiteral("(function")) ||
           body.startsWith(QStringLiteral("function"));
}

} // namespace

// =============================================================================

ScriptPropertyModel::ScriptPropertyModel(QObject* parent)
    : StructTableModel(parent) { }


// =============================================================================

ScriptConsole::ScriptConsole(QObject* p) : QObject(p) { }

namespace {

QString script_value_to_string(QJSValue const& value) {
    if (value.isUndefined()) return {};
    if (value.isNull()) return QStringLiteral("null");
    return value.toString();
}

QString format_script_log(QJSValue const& a,
                          QJSValue const& b,
                          QJSValue const& c,
                          QJSValue const& d,
                          QJSValue const& e,
                          QJSValue const& f,
                          QJSValue const& g,
                          QJSValue const& h) {
    QStringList parts;

    for (auto const& value : { a, b, c, d, e, f, g, h }) {
        if (!value.isUndefined()) parts << script_value_to_string(value);
    }

    return parts.join(' ');
}

} // namespace

void ScriptConsole::log(QJSValue a,
                        QJSValue b,
                        QJSValue c,
                        QJSValue d,
                        QJSValue e,
                        QJSValue f,
                        QJSValue g,
                        QJSValue h) {
    auto msg = format_script_log(a, b, c, d, e, f, g, h);
    qInfo().noquote() << "[Script]" << msg;
    emit logged((int)ScriptLogLevel::Log, msg);
}

void ScriptConsole::warn(QJSValue a,
                         QJSValue b,
                         QJSValue c,
                         QJSValue d,
                         QJSValue e,
                         QJSValue f,
                         QJSValue g,
                         QJSValue h) {
    auto msg = format_script_log(a, b, c, d, e, f, g, h);
    qWarning().noquote() << "[Script]" << msg;
    emit logged((int)ScriptLogLevel::Warn, msg);
}

void ScriptConsole::error(QJSValue a,
                          QJSValue b,
                          QJSValue c,
                          QJSValue d,
                          QJSValue e,
                          QJSValue f,
                          QJSValue g,
                          QJSValue h) {
    auto msg = format_script_log(a, b, c, d, e, f, g, h);
    qCritical().noquote() << "[Script]" << msg;
    emit logged((int)ScriptLogLevel::Error, msg);
}


// =============================================================================

static QStringList get_builtin_scripts() {
    QStringList files;

    QDirIterator it(":/assets/scripts/", // e.g. ":/scripts"
                    QStringList() << "*.js",
                    QDir::Files,
                    QDirIterator::Subdirectories);

    while (it.hasNext()) {
        files << it.next();
    }

    return files;
}

Script::Script(QObject* parent)
    : QObject { parent }, m_properties { new ScriptPropertyModel(this) } {
    connect(this, &Script::code_changed, this, &Script::parse);

    set_working_directory(QDir::homePath());
    set_builtin_scripts(get_builtin_scripts());
}

void Script::set_database(db::Database* db) {
    qDebug() << Q_FUNC_INFO << db;
    m_database = db;
}

bool Script::parse() {
    QStringList             errors;
    QVector<ScriptProperty> properties;
    QString                 title;
    QStringList             descriptions;

    auto header = first_comment_block(code());
    if (!header.error.isEmpty()) { errors << header.error; }

    bool saw_directive = false;
    bool saw_title     = false;

    if (errors.isEmpty()) {
        auto lines = header.text.split('\n');
        for (auto const& raw_line : std::as_const(lines)) {
            auto line = normalize_header_line(raw_line);
            if (line.isEmpty()) { continue; }

            auto keyword_end = line.indexOf(QRegularExpression("\\s"));
            auto keyword     = keyword_end < 0 ? line : line.left(keyword_end);
            auto rest =
                keyword_end < 0 ? QString {} : line.mid(keyword_end).trimmed();
            keyword = keyword.toUpper();

            if (!saw_directive) {
                saw_directive = true;
                if (keyword != "TITLE") {
                    errors << "Header must start with TITLE";
                }
            }

            if (keyword == "TITLE") {
                if (saw_title) {
                    errors << "Header contains multiple TITLE lines";
                }
                saw_title = true;
                title     = rest;
                if (title.isEmpty()) { errors << "TITLE requires text"; }
                continue;
            }

            if (keyword == "DESC") {
                descriptions << rest;
                continue;
            }

            if (keyword == "PROPERTY") {
                auto property = parse_property(line);
                if (!property.error.isEmpty()) {
                    errors << QString("PROPERTY %1: %2")
                                  .arg(property.name.isEmpty() ? "<unknown>"
                                                               : property.name,
                                       property.error);
                }
                properties << property;
                continue;
            }

            errors << QString("Unknown header directive: %1").arg(keyword);
        }
    }

    if (!saw_title && errors.isEmpty()) {
        errors << "Header must include TITLE";
    }

    set_title(title);
    set_description(descriptions.join('\n'));
    set_parse_errors(errors);
    set_valid(errors.isEmpty());
    m_properties->replace(properties);

    return valid();
}

void Script::run() {
    qDebug() << Q_FUNC_INFO;

    if (!m_database) {
        auto message = QStringLiteral("No database available.");
        qCritical().noquote() << "[Script]" << message;
        emit logged((int)ScriptLogLevel::Error, message);
        qDebug() << Q_FUNC_INFO << "No db.";
        return;
    }

    // sync for now
    auto engine = std::make_unique<QJSEngine>();
    auto api    = std::make_unique<ScriptDBInterface>(m_database);
    api->update_working_directory(working_directory());

    QStringList stack_trace;
    engine->installExtensions(QJSEngine::TranslationExtension |
                              QJSEngine::GarbageCollectionExtension);

    auto js_api_obj = engine->newQObject(api.get());
    engine->globalObject().setProperty("db", js_api_obj);

    auto console = new ScriptConsole(engine.get());
    connect(console, &ScriptConsole::logged, this, &Script::logged);
    engine->globalObject().setProperty("console", engine->newQObject(console));

    auto header = first_comment_block(code());
    auto source = looks_like_legacy_function_script(code(), header.body_start)
                      ? code()
                      : script_function(code(), header.body_start,
                                        m_properties->vector());

    auto object = engine->evaluate(source, title(), 1, &stack_trace);

    if (!stack_trace.isEmpty()) {
        qDebug() << Q_FUNC_INFO << object.toString();
        auto exception = QString("Script evaluation exception: %1")
                             .arg(object.property("name").toString());
        auto line = QString("Line number %1: %2")
                        .arg(object.property("lineNumber").toInt())
                        .arg(object.toString());
        qCritical().noquote() << "[Script]" << exception;
        qCritical().noquote() << "[Script]" << line;
        emit logged((int)ScriptLogLevel::Error, exception);
        emit logged((int)ScriptLogLevel::Error, line);
        return;
    }

    if (!object.isCallable()) {
        auto message = QStringLiteral(
            "Script did not produce a callable function.");
        qCritical().noquote() << "[Script]" << message;
        emit logged((int)ScriptLogLevel::Error, message);
        return;
    }

    // collect args

    QJSValueList list;

    for (auto const& arg : *m_properties) {
        bool ok = false;
        if (arg.type == "integer") {
            list << QJSValue(arg.value.toInt(&ok));
        } else if (arg.type == "real") {
            list << QJSValue(arg.value.toDouble(&ok));
        } else if (arg.type == "string") {
            list << QJSValue(arg.value);
        } else if (arg.type == "vec3") {
            auto value = arg.value;
            value.replace('{', ' ');
            value.replace('}', ' ');
            value.replace('[', ' ');
            value.replace(']', ' ');
            value.replace(',', ' ');

            auto parts = value.split(' ', Qt::SkipEmptyParts);

            auto dest = engine->newArray(3);

            switch (parts.size()) {
            case 1: {
                auto d = parts.value(0).toDouble(&ok);

                if (!ok) { break; }

                dest.setProperty(0, d);
                dest.setProperty(1, d);
                dest.setProperty(2, d);
                list << dest;
                break;
            }
            case 3: {
                auto d1 = parts.value(0).toDouble(&ok);
                if (!ok) { break; }
                auto d2 = parts.value(1).toDouble(&ok);
                if (!ok) { break; }
                auto d3 = parts.value(2).toDouble(&ok);
                if (!ok) { break; }

                dest.setProperty(0, d1);
                dest.setProperty(1, d2);
                dest.setProperty(2, d3);
                list << dest;
                break;
            } break;
            default: break;
            }
        }

        if (!ok) {
            qDebug() << Q_FUNC_INFO << "Bad arg";
            auto message =
                QString("Invalid argument value for %1.").arg(arg.name);
            qCritical().noquote() << "[Script]" << message;
            emit logged((int)ScriptLogLevel::Error, message);
            return;
        }
    }

    auto call_ret = object.call(list);

    if (call_ret.isError()) {
        qDebug() << Q_FUNC_INFO << "Bad call" << call_ret.toString();
        auto exception = QString("Script evaluation exception: %1")
                             .arg(call_ret.property("name").toString());
        auto line = QString("Line number %1: %2")
                        .arg(call_ret.property("lineNumber").toInt())
                        .arg(call_ret.toString());
        qCritical().noquote() << "[Script]" << exception;
        qCritical().noquote() << "[Script]" << line;
        emit logged((int)ScriptLogLevel::Error, exception);
        emit logged((int)ScriptLogLevel::Error, line);
    }
}

void Script::notify_error(QString message) {
    emit notify(ANotification::error(std::move(message)));
}

QString Script::api_markdown() {
    ScriptDBInterface api(m_database);
    return SchemaBuilder::build_markdown(&api);
}

} // namespace SolTrace::GUI::Script
