#pragma once

#include <QObject>
#include <QQuick3DGeometry>
#include <QtQuick3D/QQuick3DInstancing>

#include "job_run_common.h"
#include "job_run_thread.h"
#include "simulation_result.hpp"
#include "utilities/qt_helpers.h"

namespace SD = SolTrace::Data;
namespace RD = SolTrace::Result;

enum class RunType {
    Thread,
    Process,
};


/// Models a running simulation.
///
/// Provides pause and resume (if the simulation supports it)
/// Supports progress percent and text
///
/// When finished, users can collect results using the `take()` function.
/// When done (either finished or errored out), this object will destroy itself.
class RunningJob : public QObject {
    Q_OBJECT

    void* m_watcher;

    std::shared_ptr<db::SimulationResult> m_result;

public:
    explicit RunningJob(
        SimDataPtr          data,
        RunType             type,
        uint32_t            thread_count,
        ThreadRunnerBackend backend = ThreadRunnerBackend::Native,
        QObject*            parent  = nullptr);
    virtual ~RunningJob();

    std::shared_ptr<db::SimulationResult> take();

public slots:
    void pause();
    void resume();
    void cancel();

signals:
    void progress_update(int);
    void progress_text_update(QString);
    void finished();
    void error(QString);
};
