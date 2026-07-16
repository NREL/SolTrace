#include "schema_builder.h"
#include "script/script_db_interface.h"

#include <QDebug>
#include <QFile>
#include <QJsonArray>
#include <QJsonObject>
#include <QMetaMethod>
#include <QMetaObject>
#include <QStringList>

#include <algorithm>

namespace SolTrace::GUI::Script {

struct FunctionExport {
    QString name;
    QString desc;
};

static QString load_export_docs() {
    QFile file(QStringLiteral(":/docs/script/db_schema.md"));
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qCritical() << "Unable to load script schema docs from resource"
                    << file.fileName() << file.errorString();
        return {};
    }

    return QString::fromUtf8(file.readAll());
}

static QVector<FunctionExport> parse_exports(QString const& docs) {
    QVector<FunctionExport> ret;
    auto                    lines = docs.split('\n');
    int                     pos   = 0;

    while (pos < lines.size()) {
        while (pos < lines.size() && lines[pos].trimmed() != "---") {
            ++pos;
        }
        if (pos >= lines.size()) break;
        ++pos;

        while (pos < lines.size() && lines[pos].trimmed().isEmpty()) {
            ++pos;
        }
        if (pos >= lines.size()) break;

        auto name = lines[pos].trimmed();
        ++pos;

        QStringList body_lines;
        while (pos < lines.size() && lines[pos].trimmed() != "---") {
            body_lines << lines[pos].trimmed();
            ++pos;
        }

        ret.push_back(FunctionExport {
            name,
            body_lines.join('\n').trimmed(),
        });
    }

    return ret;
}

struct FunctionRecord {
    QString                          name;
    QString                          desc;
    QVector<QPair<QString, QString>> args;
    QString                          return_type;

    QJsonObject to_object() const {
        QJsonObject ret;
        ret[QStringLiteral("name")]        = name;
        ret[QStringLiteral("description")] = desc;

        if (args.size()) {
            QJsonArray jargs;

            for (auto const& p : args) {
                jargs.push_back(QJsonObject {
                    { "name", p.first },
                    { "type", p.second },
                });
            }

            ret[QStringLiteral("args")] = jargs;
        }

        if (return_type.size()) {
            ret[QStringLiteral("returns")] = return_type;
        }

        return ret;
    }
};

static QString clean_type(QString s) {
    s.replace("void", "");
    s.replace("db::", "");
    s.replace("QVector<Entity>", "[Entity]");
    s.replace("QList<", "[");
    s.replace(">", "]");
    return s;
}

QJsonObject SchemaBuilder::build(ScriptDBInterface* iface, QString task) {

    QJsonArray all_methods;

    QHash<QString, QMetaMethod> method_lookup;

    auto meta_obj = iface->metaObject();

    for (int i = 0; i < meta_obj->methodCount(); i++) {
        auto meta_func                  = meta_obj->method(i);
        method_lookup[meta_func.name()] = meta_func;
    }

    auto exports = parse_exports(load_export_docs());

    for (auto const& exp : exports) {
        auto name = exp.name;
        auto desc = exp.desc;
        if (!method_lookup.contains(name)) {
            qCritical() << "Missing script function" << name << "for schema";
        }

        auto const& info = method_lookup.value(name);

        FunctionRecord record;

        record.name = name;
        record.desc = desc;

        auto pnames = info.parameterNames();
        auto ptypes = info.parameterTypes();

        for (int i = 0; i < std::min(pnames.size(), ptypes.size()); ++i) {
            record.args.emplace_back(clean_type(pnames[i]),
                                     clean_type(ptypes[i]));
        }

        record.return_type = clean_type(info.typeName());

        all_methods.push_back(record.to_object());
    }

    return QJsonObject {
        { QStringLiteral("app"), QStringLiteral("SolTrace") },
        { QStringLiteral("runtime"),
          QJsonObject {
              { QStringLiteral("language"), QStringLiteral("javascript") },
              { QStringLiteral("global_objects"), QJsonArray() << "db" },
          } },

        {
            QStringLiteral("api"),
            QJsonObject {
                { QStringLiteral("db"), all_methods },
            },
        },

        { QStringLiteral("task"), task },
    };
}

QString SchemaBuilder::build_markdown(ScriptDBInterface* iface, QString task) {

    auto schema = build(iface, task);

    QStringList lines;
    lines << QStringLiteral("# SolTrace Script API") << QString {}
          << QStringLiteral("Runtime: JavaScript")
          << QStringLiteral("Global objects: db");

    if (!task.isEmpty()) {
        lines << QString {} << QStringLiteral("Task: %1").arg(task);
    }

    auto methods = schema[QStringLiteral("api")]
                       .toObject()[QStringLiteral("db")]
                       .toArray();

    for (auto const& value : methods) {
        auto method = value.toObject();

        QStringList args;
        for (auto const& arg_value : method[QStringLiteral("args")].toArray()) {
            auto arg  = arg_value.toObject();
            auto name = arg[QStringLiteral("name")].toString();
            auto type = arg[QStringLiteral("type")].toString();

            args << (type.isEmpty() ? name
                                    : QStringLiteral("%1: %2").arg(name, type));
        }

        auto returns = method[QStringLiteral("returns")].toString();
        auto suffix  = returns.isEmpty()
                           ? QString {}
                           : QStringLiteral(" -> %1").arg(returns);

        lines << QString {}
              << QStringLiteral("## db.%1(%2)%3")
                     .arg(method[QStringLiteral("name")].toString(),
                          args.join(QStringLiteral(", ")), suffix)
              << QString {} << method[QStringLiteral("description")].toString();
    }

    return lines.join('\n').trimmed();
}

} // namespace SolTrace::GUI::Script
