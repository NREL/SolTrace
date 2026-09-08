import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STFormPanel {
    id: root

    property var module
    property bool singleColumn: false
    property int labelAlignment: Qt.AlignRight | Qt.AlignVCenter

    Label {
        Layout.fillWidth: true

        visible: root.module.current_material_name.length === 0
                 || root.module.current_geometry_name.length === 0
        text: "A material and geometry must be assigned for the element to be visible."
        color: Material.color(Material.Yellow)
        wrapMode: Label.WrapAtWordBoundaryOrAnywhere
    }

    STFormRow {
        label: "Parent"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STButton {
            Layout.fillWidth: true
            property string parent_name: root.module.parent_name
            text: parent_name.length ? parent_name : "Unassigned"
            onClicked: parent_pop.open()
            SelectElementPopup {
                id: parent_pop
                exclude: [root.module.entity]
                allowNothing: true
                onSelectedElement: (element) => root.module.parent = element
                onSelectedNothing: root.module.clear_parent()
            }
        }
    }

    STFormRow {
        label: "Material"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STButton {
            Layout.fillWidth: true
            property string material_name: root.module.current_material_name
            text: material_name.length ? material_name : "Unassigned"
            onClicked: material_pop.open()
            SelectItemPopup {
                id: material_pop
                source_model: AppData.materials.materials_list
                onSelectedEntity: (entity) => root.module.current_material = entity
            }
        }
    }

    STFormRow {
        label: "Geometry"
        singleColumn: root.singleColumn
        labelAlignment: root.labelAlignment

        STButton {
            Layout.fillWidth: true
            property string geometry_name: root.module.current_geometry_name
            text: geometry_name.length ? geometry_name : "Unassigned"
            onClicked: geometry_pop.open()
            SelectItemPopup {
                id: geometry_pop
                source_model: AppData.materials.geometry_list
                onSelectedEntity: (entity) => root.module.current_geometry = entity
            }
        }
    }
}
