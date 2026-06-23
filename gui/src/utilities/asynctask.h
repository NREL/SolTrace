#pragma once

#include <QFutureWatcherBase>
#include <QObject>
#include <QtConcurrent/qtconcurrentrun.h>
#include <functional>

class AsyncTaskBase : public QObject {
    Q_OBJECT
protected:
    bool m_failed = false;

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

/// An ephemeral async task. Use this class to launch Qt Concurrent Run capable
/// functions. Results will be collected and when the task is complete,
/// delivered to a task host.
///
/// Upon completion (either success or failure), this
/// object will delete itself.
///
/// When finished, calls sink->callback(identity, QVector<T>);
///  or calls sink->callback(identity, T);
/// When failed or cancelled, calls sink->callback(identity, QString);
template <class I, class T>
class AsyncTask : public AsyncTaskBase {
    QVector<T> m_ready;

public:
    template <class U,
              class Finished,
              class Failed,
              class Function,
              class... Args>
    AsyncTask(I          identity,
              U*         sink,
              Finished&& finished_func,
              Failed&&   failed_func,
              Function&& f,
              Args&&... args)
        : AsyncTaskBase(sink) {

        auto fail_once =
            [this, sink, identity, failed_func = std::move(failed_func)](
                QString message) {
                qDebug() << Q_FUNC_INFO << "Failure handler" << message;
                if (this->m_failed) return;

                this->m_failed = true;
                std::invoke(failed_func, sink, identity, message);
                emit this->failed();
            };

        auto future = QtConcurrent::run(f, std::forward<Args>(args)...)
                          .onFailed(sink,
                                    [fail_once](std::exception const& e) {
                                        fail_once(QString(e.what()));
                                        return T();
                                    })
                          .onFailed(sink,
                                    [fail_once](QException const& e) {
                                        fail_once(QString(e.what()));
                                        return T();
                                    })
                          .onFailed(sink, [fail_once] {
                              fail_once("unknown exception");
                              return T();
                          });

        using FW = QFutureWatcher<T>;

        auto watcher = new FW(this);

        // We use self-owning tasks
        connect(watcher, &FW::finished, this, &AsyncTask::deleteLater);

        connect(watcher,
                &FW::finished,
                sink,
                [this,
                 watcher,
                 finished_func = std::move(finished_func),
                 sink,
                 identity]() {
                    qDebug() << Q_FUNC_INFO << watcher << this;

                    if (watcher->isCanceled() or this->m_failed) { return; }

                    if constexpr (std::invocable<decltype(finished_func),
                                                 U*,
                                                 I,
                                                 QVector<T>>) {

                        std::invoke(
                            finished_func, sink, identity, std::move(m_ready));

                    } else if constexpr (std::invocable<decltype(finished_func),
                                                        U*,
                                                        I,
                                                        T>) {

                        if (m_ready.size()) {
                            std::invoke(finished_func,
                                        sink,
                                        identity,
                                        std::move(m_ready[0]));
                        } else {
                            std::invoke(finished_func, sink, identity, T());
                        }


                    } else {
                        []<bool flag = false>() {
                            static_assert(flag, "no match");
                        }();
                    }

                    qDebug() << Q_FUNC_INFO << "Finished";

                    emit this->finished();
                });

        connect(watcher, &FW::canceled, sink, [fail_once]() {
            fail_once("cancelled");
        });


        connect(watcher, &FW::resultReadyAt, this, [this, watcher](int index) {
            auto result = watcher->resultAt(index);

            m_ready << result;
        });

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

    bool has_results() const { return !m_ready.empty(); }

    /// Obtain all results collected at this point
    QVector<T> all_results() const { return m_ready; }

    /// Obtain the first result collected. If no results have been collected
    /// (check has_results), returns a default constructed value.
    T result() const { return m_ready.value(0); }
};


/// Launch a task in another thread.
///
/// Most standalone functions are supported; special support remains for
/// QPromise as the first argument.
///
/// When finished, calls sink->callback(identity, QVector<T>);
///  or calls sink->callback(identity, T);
/// When failed or cancelled, calls sink->callback(identity, QString);
///
/// NOTE! for proper handling, please make sure to always emplace a result.
template <class T,
          class I,
          class U,
          class Finished,
          class Failed,
          class Function,
          class... Args>
AsyncTask<I, T>* launch_async_task(I          identity,
                                   U*         sink,
                                   Finished&& finished_func,
                                   Failed&&   failed_func,
                                   Function&& f,
                                   Args&&... args) {
    return new AsyncTask<I, T>(identity,
                               sink,
                               std::forward<Finished>(finished_func),
                               std::forward<Failed>(failed_func),
                               f,
                               std::forward<Args>(args)...);
}
