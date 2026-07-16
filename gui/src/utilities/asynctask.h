#pragma once

#include "utilities/result.h"

#include <QFutureWatcherBase>
#include <QObject>
#include <QtConcurrent/qtconcurrentrun.h>

#include <functional>
#include <type_traits>

class AsyncTaskBase : public QObject {
    Q_OBJECT

signals:
    void internal_cancel(QPrivateSignal);

public:
    explicit AsyncTaskBase(QObject* parent);
    ~AsyncTaskBase() override;

public slots:
    void cancel();

signals:
    /// Successfully completed.
    void finished();

    /// Was not able to complete or was cancelled
    void failed();

    /// 0 - 100
    void progress(int);

    void progress_status(QString);
};

struct TaskControl {
    TaskControl();

    TaskControl(TaskControl const&) = delete;
    TaskControl(TaskControl&&)      = delete;

    TaskControl& operator=(TaskControl const&) = delete;
    TaskControl& operator=(TaskControl&&)      = delete;

    virtual void suspendIfRequested() const = 0;
    virtual bool cancelRequested() const    = 0;

    virtual void setProgressValue(int) const                 = 0;
    virtual void setProgressValueAndText(int, QString) const = 0;
};

template <class SuccessType, class FailureType>
struct TaskControlDerived : public TaskControl {
    QPromise<Result<SuccessType, FailureType>>& promise;

    TaskControlDerived(QPromise<Result<SuccessType, FailureType>>& _promise)
        : promise(_promise) { }

    void suspendIfRequested() const override {
        return promise.suspendIfRequested();
    }
    bool cancelRequested() const override { return promise.isCanceled(); }

    virtual void setProgressValue(int value) const override {
        promise.setProgressValue(value);
    }
    virtual void setProgressValueAndText(int     value,
                                         QString text) const override {
        promise.setProgressValueAndText(value, text);
    }
};

/// An ephemeral async task. Use this class to launch Qt Concurrent Run capable
/// functions. Results will be collected and when the task is complete,
/// delivered to a task host.
///
/// Upon completion (either success or failure), this
/// object will delete itself.
///
/// When finished, calls sink->callback(identity, SuccessType);
/// When failed or cancelled, calls sink->callback(identity, FailureType);
///
/// The failure type should have a conversion from a QString to record
/// exceptions.
template <class IdentityType, class SuccessType, class FailureType>
class AsyncTask : public AsyncTaskBase {
    // QVector<Result<SuccessType, FailureType>> m_ready;

public:
    template <class SinkType,
              class FinishedFunc,
              class FailedFunc,
              class TaskFunction,
              class... Args>
    AsyncTask(IdentityType   identity,
              SinkType*      sink,
              FinishedFunc&& finished_func,
              FailedFunc&&   failed_func,
              TaskFunction&& f,
              Args&&... args)
        : AsyncTaskBase(sink) {

        using ResultType = Result<SuccessType, FailureType>;

        auto task_function = [f, identity](QPromise<ResultType>& promise,
                                           std::decay_t<Args>... args) {
            try {
                auto control =
                    TaskControlDerived<SuccessType, FailureType> { promise };

                // hide details...
                TaskControl& dispatch = control;

                promise.setProgressRange(0, 100);
                ResultType result = f(dispatch, std::move(args)...);
                promise.emplaceResult(std::move(result));
                qDebug() << Q_FUNC_INFO << identity << "Success";

            } catch (std::exception const& e) {
                qDebug() << Q_FUNC_INFO << "Failure handler" << e.what();
                promise.emplaceResult(FailureType(QString::fromUtf8(e.what())));
            } catch (...) {
                qDebug() << Q_FUNC_INFO << "Failure handler: unknown exception";
                promise.emplaceResult(
                    FailureType(QStringLiteral("unknown exception")));
            }
        };

        auto future =
            QtConcurrent::run(task_function, std::forward<Args>(args)...);

        using FW = QFutureWatcher<ResultType>;

        auto watcher = new FW(this);

        // We use self-owning tasks
        connect(watcher, &FW::finished, this, &AsyncTask::deleteLater);

        auto on_completion = [this,
                              watcher,
                              finished_func = std::move(finished_func),
                              failed_func   = std::move(failed_func),
                              sink,
                              identity]() {
            qDebug() << Q_FUNC_INFO << watcher << this;

            if (watcher->isCanceled()) {
                std::invoke(failed_func,
                            sink,
                            identity,
                            FailureType(QStringLiteral("Canceled")));
                emit this->failed();
                return;
            }

            // we _should_ have a result here, and only one
            auto result = watcher->future().takeResult();

            if (result.is_failure()) {
                std::invoke(failed_func,
                            sink,
                            identity,
                            std::move(result.get_failure()));
                emit this->failed();
                return;
            }

            // we _should_ be success here

            std::invoke(
                finished_func, sink, identity, std::move(result.get_success()));

            qDebug() << Q_FUNC_INFO << "Finished";

            emit this->finished();
        };

        connect(watcher, &FW::finished, sink, on_completion);

        // Permit progress watching
        connect(watcher, &FW::progressValueChanged, this, [this](int value) {
            emit progress(value);
        });

        connect(watcher, &FW::progressTextChanged, this, [this](QString value) {
            emit progress_status(value);
        });

        // Set up cancelling
        connect(this, &AsyncTask::internal_cancel, watcher, &FW::cancel);

        watcher->setFuture(future);
    }

    ~AsyncTask() override { }
};


/// Launch a task in another thread.
///
/// When finished, calls sink->callback(identity, SuccessType);
/// When failed or cancelled, calls sink->callback(identity, FailureType);
///
/// The failure type should have a conversion from a QString to record
/// exceptions.
///
/// By default the first argument must be TaskControl &, which you can use to
/// check for cancellation or suspension.
template <class SuccessType,
          class FailureType,
          class IdentityType,
          class SinkType,
          class FinishedFunc,
          class FailedFunc,
          class TaskFunction,
          class... Args>
AsyncTask<IdentityType, SuccessType, FailureType>*
launch_async_task(IdentityType   identity,
                  SinkType*      sink,
                  FinishedFunc&& finished_func,
                  FailedFunc&&   failed_func,
                  TaskFunction&& f,
                  Args&&... args) {
    return new AsyncTask<IdentityType, SuccessType, FailureType>(
        identity,
        sink,
        std::forward<FinishedFunc>(finished_func),
        std::forward<FailedFunc>(failed_func),
        f,
        std::forward<Args>(args)...);
}

/// Helper macro
#define ASYNC_TASK_SYNC_POINT(CONTROL)                                         \
    {                                                                          \
        CONTROL.suspendIfRequested();                                          \
        if (CONTROL.cancelRequested()) {                                       \
            return return_failure(QStringLiteral("Cancelled"));                \
        }                                                                      \
    }
