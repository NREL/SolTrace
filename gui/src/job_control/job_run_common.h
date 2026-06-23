#pragma once

#include <QString>
#include <QtCore/qpromise.h>

#include "database/database.h"
#include "database/simulationresult.h"

#include "simulation_data_api.hpp"

namespace SD     = SolTrace::Data;
using SimDataPtr = std::shared_ptr<db::DatabaseExport>;


using ResultPtr = std::shared_ptr<db::SimulationResult>;
using SimResult = std::variant<ResultPtr, QString>;
