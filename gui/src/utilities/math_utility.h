#pragma once

#include "simulation_data_api.hpp"

#include <glm/gtc/quaternion.hpp>

#include <QDebug>
#include <QQuaternion>
#include <QVector3D>

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

/*!
 * \brief Compute linear interpolation
 * \param x The interpolation source.
 * \param x0 The minimum range of x.
 * \param x1 The maximum range of x.
 * \param y0 The minimum range of y.
 * \param y1 The maximum range of y.
 *
 * Template T and U should support +, -, *, and / operators.
 */
template <class T, class U>
U lerp(T x, T const& x0, T const& x1, U const& y0, U const& y1) {
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
}


// NOTE that we are truncating vectors to floats...
namespace SD = SolTrace::Data;

inline QVector3D convert(glm::dvec3 v) {
    return QVector3D(v.x, v.y, v.z);
}

inline glm::dvec3 convert(QVector3D v) {
    return { v.x(), v.y(), v.z() };
}


inline QQuaternion convert(glm::dquat v) {
    return QQuaternion(v.w, v.x, v.y, v.z);
}

inline glm::dquat convert(QQuaternion v) {
    return { v.scalar(), v.x(), v.y(), v.z() };
}


inline QDebug operator<<(QDebug debug, glm::dvec3 const& type) {
    QDebugStateSaver saver(debug);
    debug.nospace() << "dvec3(" << type.x << ", " << type.y << ", " << type.z
                    << ")";
    return debug;
}

inline QDebug operator<<(QDebug debug, glm::dquat const& type) {
    QDebugStateSaver saver(debug);
    debug.nospace() << "dquat(" << type.x << ", " << type.y << ", " << type.z
                    << ", " << type.w << ")";
    return debug;
}
