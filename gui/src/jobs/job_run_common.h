#pragma once

#include <QString>
#include <QtCore/qpromise.h>

#include "data/database.h"
#include "data/simulationresult.h"

#include "simulation_data_api.hpp"

namespace SD = SolTrace::Data;

namespace SolTrace::GUI::Jobs {

/// Exported simulation input and source database retained for result mapping.
using SimDataPtr = std::shared_ptr<SolTrace::GUI::Data::DatabaseExport>;

/// GUI-owned simulation result pointer.
using ResultPtr = std::shared_ptr<SolTrace::GUI::Data::SimulationResult>;

/// Thread runner output: a result pointer on success or an error message.
using SimResult = std::variant<ResultPtr, QString>;

} // namespace SolTrace::GUI::Jobs
