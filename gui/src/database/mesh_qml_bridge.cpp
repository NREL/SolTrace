#include "mesh_qml_bridge.h"

#include "database/geometryeditor.h"

#include <glm/common.hpp>

namespace db {

void QMLMesh::rebuild_geometry() {
    // Always clear at start, remove stale geom
    clear();

    // If mesh is empty, bail
    if (m_current_mesh.vertex.empty() || m_current_mesh.triangles.empty()) {
        qWarning() << Q_FUNC_INFO
                   << "Geometry is empty, or unable to be generated";
        update(); // Signal empty
        return;
    }

    // Compute AABB
    constexpr float max_float = std::numeric_limits<float>::max();

    glm::vec3 bounds_min(max_float);
    glm::vec3 bounds_max(-max_float);

    for (auto const& p : std::as_const(m_current_mesh.vertex)) {
        bounds_min = glm::min(bounds_min, p.position);
        bounds_max = glm::max(bounds_max, p.position);
    }

    // Init fresh buffers
    auto indexBuffer = QByteArray(
        reinterpret_cast<const char*>(m_current_mesh.triangles.data()),
        m_current_mesh.triangles.size() * sizeof(glm::uvec3));
    auto vertexBuffer =
        QByteArray(reinterpret_cast<const char*>(m_current_mesh.vertex.data()),
                   m_current_mesh.vertex.size() * sizeof(Vertex));

    addAttribute(QQuick3DGeometry::Attribute::PositionSemantic,
                 offsetof(Vertex, position),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::NormalSemantic,
                 offsetof(Vertex, normal),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::TexCoord0Semantic,
                 offsetof(Vertex, uv),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::IndexSemantic,
                 0,
                 QQuick3DGeometry::Attribute::ComponentType::U32Type);

    auto bb = BoundingBox {
        .min = QVector3D(bounds_min.x, bounds_min.y, bounds_min.z),
        .max = QVector3D(bounds_max.x, bounds_max.y, bounds_max.z),
    };

    setStride(sizeof(Vertex));
    setVertexData(vertexBuffer);
    setIndexData(indexBuffer);
    setBounds(bb.min, bb.max);
    setPrimitiveType(QQuick3DGeometry::PrimitiveType::Triangles);

    qDebug() << Q_FUNC_INFO << m_current_mesh.triangles.size()
             << m_current_mesh.vertex.size();

    update();
}

QMLMesh::QMLMesh() {
    update();

    connect(
        this, &QMLMesh::current_mesh_changed, this, &QMLMesh::rebuild_geometry);
}

QMLMesh::~QMLMesh() { }

const Mesh& QMLMesh::current_mesh() const {
    return m_current_mesh;
}

void QMLMesh::set_current_mesh(Mesh const& new_value) {
    // Always assume it is different, as an equality check is too expensive
    m_current_mesh = new_value;
    emit current_mesh_changed();
}

} // namespace db