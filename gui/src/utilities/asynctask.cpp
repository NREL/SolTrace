#include "asynctask.h"


AsyncTaskBase::AsyncTaskBase(QObject* parent) : QObject { parent } { }

AsyncTaskBase::~AsyncTaskBase() = default;

void AsyncTaskBase::cancel() {
    emit internal_cancel(QPrivateSignal {});
}
