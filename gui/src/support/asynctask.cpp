#include "asynctask.h"

namespace SolTrace::GUI::Support {

// anchor vtable
TaskControl::TaskControl() = default;

AsyncTaskBase::AsyncTaskBase(QObject* parent) : QObject { parent } { }

AsyncTaskBase::~AsyncTaskBase() = default;

void AsyncTaskBase::cancel() {
    emit internal_cancel(QPrivateSignal {});
}

} // namespace SolTrace::GUI::Support