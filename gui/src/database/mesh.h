#pragma once

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

#include <QVector>

namespace db {

/// Vertex format used by QQuick3DGeometry adapters.
struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv; // quick geometry only supports floating
};

/// Triangle mesh exchanged between GUI analysis code and Quick3D adapters.
struct Mesh {
    // using QVector for COW
    QVector<Vertex>     vertex;
    QVector<glm::uvec3> triangles;
};


} // namespace db
