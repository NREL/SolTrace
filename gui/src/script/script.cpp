#include "script.h"
#include "script/schema_builder.h"

#include <QDebug>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QRegularExpression>
#include <QUrl>

#include "modules/database_module.h"
#include "modules/simulation_module.h"

#include <initializer_list>
#include <optional>
#include <utility>

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

QString normalize_config_key(QString key) {
    return key.trimmed().toLower().replace('-', "_").replace(' ', "_");
}

QJsonValue config_value(QJsonObject const& config,
                        std::initializer_list<QString> keys) {
    for (auto const& key : keys) {
        auto value = config.value(key);
        if (!value.isUndefined()) return value;
    }

    for (auto it = config.constBegin(); it != config.constEnd(); ++it) {
        auto const input_key = normalize_config_key(it.key());
        for (auto const& key : keys) {
            if (input_key == normalize_config_key(key)) return it.value();
        }
    }

    return {};
}

bool json_bool(QJsonValue const& value, bool fallback) {
    return value.isUndefined() ? fallback : value.toBool(fallback);
}

uint32_t json_uint(QJsonValue const& value, uint32_t fallback) {
    if (value.isUndefined()) return fallback;
    auto number = value.toDouble(fallback);
    if (number < 0.0) return fallback;
    return static_cast<uint32_t>(number);
}

std::optional<App::SimulationModule::Runner>
runner_from_value(QJsonValue const& value,
                  App::SimulationRunnerModel const* runners) {
    if (value.isUndefined()) return std::nullopt;

    if (value.isDouble()) {
        return static_cast<App::SimulationModule::Runner>(value.toInt());
    }

    auto key = normalize_config_key(value.toString());
    if (key == QStringLiteral("cpu") || key == QStringLiteral("legacy") ||
        key == QStringLiteral("native")) {
        return App::SimulationModule::CPU;
    }
    if (key == QStringLiteral("embree")) return App::SimulationModule::Embree;
    if (key == QStringLiteral("gpu") || key == QStringLiteral("optix")) {
        return App::SimulationModule::GPU;
    }

    if (key.endsWith(QStringLiteral("_runner"))) {
        key.chop(QStringLiteral("_runner").size());
    }

    if (runners) {
        for (int i = 0; i < runners->rowCount(); ++i) {
            auto const* record = runners->get_at(i);
            if (!record) continue;
            auto name = normalize_config_key(record->name);
            if (name == key) return record->runner;
        }
    }

    return std::nullopt;
}

QString resolve_script_file_path(QString const& working_directory,
                                 QString const& relative_path,
                                 bool           must_exist) {
    if (relative_path.isEmpty() || QFileInfo(relative_path).isAbsolute()) {
        return {};
    }

    auto base_dir = QDir(working_directory.isEmpty() ? QDir::currentPath()
                                                     : working_directory);
    auto base_path = QDir::cleanPath(base_dir.absolutePath());
    auto file_path = QDir::cleanPath(base_dir.filePath(relative_path));

    if (file_path == base_path ||
        !file_path.startsWith(base_path + QDir::separator())) {
        return {};
    }

    if (must_exist && !QFileInfo(file_path).isFile()) return {};
    return file_path;
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

ScriptSimulationInterface::ScriptSimulationInterface(App::SimulationModule* sim,
                                                     QObject* parent)
    : QObject(parent), m_simulation(sim) { }

bool ScriptSimulationInterface::start(QJsonObject config) {
    if (!m_simulation || m_simulation->is_running()) return false;
    if (!m_simulation->current_database()) return false;

    auto old_runner         = m_simulation->runner();
    auto old_ray_count      = m_simulation->ray_count();
    auto old_max_ray_count  = m_simulation->max_ray_count();
    auto old_max_threads    = m_simulation->max_threads();
    auto old_seed_value     = m_simulation->seed_value();
    auto old_sun_shape      = m_simulation->sun_shape();
    auto old_optical_errors = m_simulation->optical_errors();
    auto old_point_focus    = m_simulation->point_focus_system();

    if (auto runner = runner_from_value(
            config_value(config,
                         { QStringLiteral("runner"),
                           QStringLiteral("runner_name") }),
            m_simulation->runners())) {
        m_simulation->set_runner(*runner);
    }

    if (auto value = config_value(config, { QStringLiteral("runner_index") });
        !value.isUndefined()) {
        m_simulation->set_runner(
            m_simulation->runners()->runner_at(value.toInt()));
    }

    m_simulation->set_ray_count(json_uint(
        config_value(config, { QStringLiteral("ray_count"),
                               QStringLiteral("rays") }),
        m_simulation->ray_count()));
    m_simulation->set_max_ray_count(json_uint(
        config_value(config, { QStringLiteral("max_ray_count"),
                               QStringLiteral("max_rays_traced") }),
        m_simulation->max_ray_count()));
    m_simulation->set_max_threads(json_uint(
        config_value(config, { QStringLiteral("max_threads"),
                               QStringLiteral("threads") }),
        m_simulation->max_threads()));
    m_simulation->set_seed_value(json_uint(
        config_value(config, { QStringLiteral("seed_value"),
                               QStringLiteral("seed") }),
        m_simulation->seed_value()));
    m_simulation->set_sun_shape(json_bool(
        config_value(config, { QStringLiteral("sun_shape"),
                               QStringLiteral("include_sun_shape_errors") }),
        m_simulation->sun_shape()));
    m_simulation->set_optical_errors(json_bool(
        config_value(config, { QStringLiteral("optical_errors"),
                               QStringLiteral("include_optical_errors") }),
        m_simulation->optical_errors()));
    m_simulation->set_point_focus_system(json_bool(
        config_value(config, { QStringLiteral("point_focus_system") }),
        m_simulation->point_focus_system()));

    m_simulation->run();
    auto started = m_simulation->is_running();

    m_simulation->set_runner(old_runner);
    m_simulation->set_ray_count(old_ray_count);
    m_simulation->set_max_ray_count(old_max_ray_count);
    m_simulation->set_max_threads(old_max_threads);
    m_simulation->set_seed_value(old_seed_value);
    m_simulation->set_sun_shape(old_sun_shape);
    m_simulation->set_optical_errors(old_optical_errors);
    m_simulation->set_point_focus_system(old_point_focus);

    return started;
}

ScriptSceneInterface::ScriptSceneInterface(App::DatabaseModule* databases,
                                           ScriptDBInterface* database_interface,
                                           QString working_directory,
                                           QObject* parent)
    : QObject(parent),
      m_databases(databases),
      m_database_interface(database_interface),
      m_working_directory(std::move(working_directory)) { }

bool ScriptSceneInterface::set_current_index(int index) {
    if (!m_databases || !m_databases->set_current(index)) return false;
    if (m_database_interface) {
        m_database_interface->set_database(m_databases->current_database());
    }
    return true;
}

bool ScriptSceneInterface::set_current_name(QString name) {
    if (!m_databases) return false;

    for (int i = 0; i < m_databases->rowCount(); ++i) {
        auto const* record = m_databases->get_at(i);
        if (record && record->database &&
            record->database->name() == name) {
            return set_current_index(i);
        }
    }

    return false;
}

bool ScriptSceneInterface::new_blank(QString name) {
    if (!m_databases) return false;
    m_databases->append_new(name.isEmpty() ? QStringLiteral("Untitled") : name);
    if (m_database_interface) {
        m_database_interface->set_database(m_databases->current_database());
    }
    return true;
}

bool ScriptSceneInterface::new_from_file(QString relative_path,
                                         QString name_override) {
    if (!m_databases) return false;

    auto path = resolve_script_file_path(m_working_directory,
                                         relative_path,
                                         true /* must_exist */);
    if (path.isEmpty()) return false;

    m_databases->load_url(QUrl::fromLocalFile(path), name_override);
    return true;
}

bool ScriptSceneInterface::export_json(QString relative_path) {
    if (!m_databases) return false;

    auto path = resolve_script_file_path(m_working_directory,
                                         relative_path,
                                         false /* must_exist */);
    if (path.isEmpty()) return false;

    return m_databases->export_current_json(path);
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

void Script::set_services(App::DatabaseModule* databases,
                          App::SimulationModule* simulation) {
    m_databases  = databases;
    m_simulation = simulation;
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
        emit runCompleted();
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

    auto sim_api = std::make_unique<ScriptSimulationInterface>(m_simulation);
    engine->globalObject().setProperty("sim",
                                       engine->newQObject(sim_api.get()));

    auto scene_api = std::make_unique<ScriptSceneInterface>(
        m_databases, api.get(), working_directory());
    engine->globalObject().setProperty("scenes",
                                       engine->newQObject(scene_api.get()));

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
        emit runCompleted();
        return;
    }

    if (!object.isCallable()) {
        auto message = QStringLiteral(
            "Script did not produce a callable function.");
        qCritical().noquote() << "[Script]" << message;
        emit logged((int)ScriptLogLevel::Error, message);
        emit runCompleted();
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
            emit runCompleted();
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

    emit runCompleted();
}

void Script::notify_error(QString message) {
    emit notify(ANotification::error(std::move(message)));
}

QString Script::api_markdown() {
    ScriptDBInterface api(m_database);
    return SchemaBuilder::build_markdown(&api);
}

} // namespace SolTrace::GUI::Script
