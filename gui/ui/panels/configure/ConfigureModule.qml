import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    Binding {
        target: bar
        property: "currentIndex"
        value: App.view.configure_section
    }

    Label {
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8

        text: "Configure the ray source, materials, and geometry. Then bind the material and geometry to an element, and stage the elements in the scene."
        wrapMode: Text.WordWrap
    }

    STPipelineBar {
        id: bar
        currentIndex: App.view.configure_section
        onCurrentIndexChanged: App.view.configure_section = currentIndex

        Layout.fillWidth: true

        collapseLabels: App.view.left_panel.size === SplitPanelData.Small

        prefixModel: ["2a", "2b", "2c", "2d"]
        iconModel: ["\uf185", "\uf042", "\uf1b2", "\uf5ee"]
        model: ["Ray Source", "Materials", "Geometries", "Staging"]
        
    }
    
    StackLayout {
        currentIndex: App.view.configure_section

        ConfigureSun {
        }

        ConfigureMaterials {
        }
        
        ConfigureGeometry {
        }

        ConfigureLayout {
        }
    }

    Rectangle {
        Layout.fillWidth: true
        height: 1
        color: Material.dividerColor
    }
    
    WorkflowStepper {
        Layout.fillWidth: true
        previous: "Load Scene"
        next: "Run Tracer"
        currentIndex: ViewModule.Configure
    }
}
