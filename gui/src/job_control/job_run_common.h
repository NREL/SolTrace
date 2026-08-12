#pragma once

#include <QString>
#include <QtCore/qpromise.h>

#include "database/database.h"
#include "database/simulationresult.h"

#include "simulation_data_api.hpp"

namespace SD = SolTrace::Data;
/// Exported simulation input and source database retained for result mapping.
using SimDataPtr = std::shared_ptr<db::DatabaseExport>;

/// GUI-owned simulation result pointer.
using ResultPtr = std::shared_ptr<db::SimulationResult>;

/// Thread runner output: a result pointer on success or an error message.
using SimResult = std::variant<ResultPtr, QString>;
