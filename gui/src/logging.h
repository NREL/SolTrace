#pragma once

#include <QObject>

#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

namespace SolTrace::GUI::App {

struct LogRecord {
    QString content;

    RECORD_META(LogRecord, SM_EXPOSE_RO(content));
};

class LogList : public StructTableModel<LogRecord> {
    Q_OBJECT

    Q_WRITABLE_PROPERTY(quint32, max_line_count, 1000);
    Q_READONLY_PROPERTY(QString, log_directory);

private slots:
    void spin_off();

public:
    LogList();
    virtual ~LogList();

public slots:
    void append_line(QString);

public:
    Q_INVOKABLE bool open_log_directory();
};

LogList* initialize_logging_handler();
LogList* current_log_list();

} // namespace SolTrace::GUI::App
