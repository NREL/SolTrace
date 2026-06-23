#pragma once

#include <QJsonObject>
#include <QString>

namespace SolTrace::GUI::Script {

class ScriptDBInterface;

class SchemaBuilder {
public:
    static QJsonObject build(ScriptDBInterface*, QString task);
};

} // namespace SolTrace::GUI::Script
