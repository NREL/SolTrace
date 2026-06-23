#pragma once

#include "database/mesh.h"

#include <QObject>
#include <QQuick3DGeometry>

namespace db {

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

    Mesh const& current_mesh() const;
    void        set_current_mesh(Mesh const& new_value);

signals:
    void current_mesh_changed();
};

} // namespace db