#pragma once

#include "database/mesh.h"

#include <QObject>
#include <QQuick3DGeometry>

namespace db {

/// QQuick3DGeometry wrapper around a db::Mesh value.
class QMLMesh : public QQuick3DGeometry {
    Q_OBJECT

    // have to be careful here. no helpers.
    Q_PROPERTY(Mesh current_mesh READ current_mesh NOTIFY current_mesh_changed)

    Mesh m_current_mesh = Mesh();

private slots:
    void rebuild_geometry();

public:
    QMLMesh();
    ~QMLMesh();

    /// Current mesh value used to rebuild the Quick3D geometry buffers.
    Mesh const& current_mesh() const;

    /// Replace current_mesh and rebuild the geometry buffers.
    void set_current_mesh(Mesh const& new_value);

signals:
    void current_mesh_changed();
};

} // namespace db
