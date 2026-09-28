#pragma once

#include <QJsonObject>
#include <QString>

namespace SolTrace::GUI::Scripting {

class ScriptDBInterface;

class SchemaBuilder {
public:
    static QJsonObject build(ScriptDBInterface*, QString task);
    static QString     build_markdown(ScriptDBInterface*, QString task = {});
};

} // namespace SolTrace::GUI::Scripting
