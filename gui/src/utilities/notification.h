#pragma once

#include <QObject>
#include <utility>

struct ANotification {
    Q_GADGET
    Q_PROPERTY(QString message MEMBER message)
    Q_PROPERTY(Type type MEMBER type)

public:
    enum Type { INFO, WARNING, ERROR };

    Q_ENUM(Type);


    QString message;
    Type    type;

    ANotification() = default;
    ANotification(QString msg, Type t = INFO)
        : message(std::move(msg)), type(t) { }

    static ANotification info(QString msg) {
        return ANotification(std::move(msg), INFO);
    }

    static ANotification warning(QString msg) {
        return ANotification(std::move(msg), WARNING);
    }

    static ANotification error(QString msg) {
        return ANotification(std::move(msg), ERROR);
    }
};

Q_DECLARE_METATYPE(ANotification)
