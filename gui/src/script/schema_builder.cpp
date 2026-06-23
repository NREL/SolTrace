#include "schema_builder.h"
#include "script/script_db_interface.h"

#include <QJsonArray>
#include <QJsonObject>
#include <QMetaMethod>
#include <QMetaObject>


namespace SolTrace::GUI::Script {

struct FunctionExport {
    const char* name;
    const char* desc;
};

static constexpr FunctionExport exports[] = {
    {
        "create",
        "Create a new CSP element. Requires geometry and material to be "
        "useful.",
    },
    {
        "get_all_materials",
        "Obtain all available materials, identified by internal IDs.",
    },
};

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
    s.replace("QList<", "[");
    s.replace(">", "]");
    return s;
}

template <class U, class F>
QVector<U> map(QVector<QByteArray> in, F&& func) {
    QVector<U> ret;

    for (auto const& item : in) {
        ret << func(item);
    }
    return ret;
}

QJsonObject SchemaBuilder::build(ScriptDBInterface* iface, QString task) {

    QJsonArray all_methods;

    QHash<QString, QMetaMethod> method_lookup;

    auto meta_obj = iface->metaObject();

    for (int i = 0; i < meta_obj->methodCount(); i++) {
        auto meta_func                  = meta_obj->method(i);
        method_lookup[meta_func.name()] = meta_func;
    }

    for (auto exp : exports) {
        auto name = QString(exp.name);
        auto desc = QString(exp.desc);
        if (!method_lookup.contains(name)) {
            qCritical() << "Missing script function" << name << "for schema";
        }

        auto const& info = method_lookup.value(name);

        FunctionRecord record;

        record.name = name;

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
              { QStringLiteral("entrypoint"),
                QStringLiteral("script evaluates to a function") },
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

} // namespace SolTrace::GUI::Script
