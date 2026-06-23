#pragma once

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

#include <QVector>

namespace db {

struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv; // quick geometry only supports floating
};

struct Mesh {
    // using QVector for COW
    QVector<Vertex>     vertex;
    QVector<glm::uvec3> triangles;
};


} // namespace db