import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    readonly property var backendFeatures: [
        "gpu_ray_tracing",
        "large_scale_simulations",
        "open_source_extensible",
        "sam_integration",
        "python_api"
    ]
    readonly property var guiFeatures: [
        "viewport",
        "geometry_editor",
        "flux_visualization",
        "theming_accessibility",
        "documentation"
    ]

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        Header {
            text: App.docs.get("settings.features", "title")
        }

        FeatureSection {
            title: App.docs.get("settings.features.backend", "title")
            keys: root.backendFeatures
            Layout.fillWidth: true
        }

        FeatureSection {
            title: App.docs.get("settings.features.gui", "title")
            keys: root.guiFeatures
            Layout.fillWidth: true
        }

        Item { Layout.fillHeight: true }
    }
}
