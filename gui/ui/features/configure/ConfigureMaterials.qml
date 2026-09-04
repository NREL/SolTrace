import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

AdaptiveFilteredEditor {
    id: root

    source_model: App.materials.materials_list

    wideThreshold: 500
    listWidth: 250
    emptyListText: qsTr("Materials define optical properties. Create one below to assign it to elements.")

    onEditingChanged: {
        App.view.editing_material = editing
    }

    Connections {
        target: App.view
        function onEditing_material_changed() {
            root.editing = App.view.editing_material
        }
    }

    Connections {
        target: App.materials
        function onCurrent_material_changed() {
            root.selectSourceIndex(
                        root.source_model.index_of(App.materials.current_material))
        }
    }

    onItemClicked: function(index, modelData) {
        // if (App.db) App.db.clear_selection()
        App.materials.current_material = modelData.entity
        // if (App.db) App.db.select_all_with_material(modelData.entity)
    }

    listHeader: RowLayout {
        STSearchField {
            id: search_field
            Layout.fillWidth: true

            Binding {
                target: root
                property: "filterText"
                value: search_field.text
            }
        }

        STIconButton {
            icon: "\uf05b"
            toolTip: qsTr("Select Material From View")
            onClicked: App.view.mouse_mode = ViewModule.SelectMaterial
        }
    }

    listFooter: RowLayout {
        CreateNewItemButton {
            title: qsTr("New Material")

            onCreateRequested: function(name) {
                var new_name = AppData.current_database.sanitize_material_name(name)
                AppData.current_database.add_material_group(new_name)
            }
        }
    }

    listDelegate: STItemDelegate {
        text: itemModel ? itemModel.name : qsTr("No name")
        highlighted: isCurrent
    }

    detailView: ColumnLayout {
        spacing: 8

        RowLayout {
            STIconButton {
                icon: "\uf053"
                visible: !root.wideMode
                onClicked: root.goBack()
            }

            RenameLabel {
                text: App.materials.current_material_name

                onAccepted: (new_name) => {
                    var mats = App.materials;
                    var curr_db = mats.current_database

                    new_name = curr_db.sanitize_material_name(new_name);

                    curr_db.set_name_of(
                                    mats.current_material,
                                    new_name
                                    )
                }
            }
        }

        ScrollView {
            id: mat_scroll
            Layout.fillHeight: true
            Layout.fillWidth: true
            contentWidth: availableWidth

            ColumnLayout {
                width: mat_scroll.availableWidth

                InlineDocumentation {
                    key: "configure.materials"
                }

                SurfacePropertyGraphic {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 148

                    property var face: bar.currentIndex == 0
                                       ? App.materials.material_edit.front_editor
                                       : App.materials.material_edit.back_editor

                    reflectance: face.reflectivity
                    transmittance: face.transmissivity
                    nFront: face.refraction_index_front
                    nBack: face.refraction_index_back
                    slopeErrorMrad: face.slope_error
                    specularityErrorMrad: face.specularity_error
                }

                STComboBar {
                    id: bar
                    Layout.fillWidth: true
                    model: [qsTr("Front"), qsTr("Back")]
                }

                SwipeView {
                    Layout.fillWidth: true
                    Layout.preferredHeight: currentItem ? currentItem.implicitHeight : 0
                    interactive: false
                    clip: true
                    currentIndex: bar.currentIndex

                    MaterialOpticals {
                        collapsed: false
                        Layout.fillWidth: true
                        side_editor: App.materials.material_edit.front_editor
                    }

                    MaterialOpticals {
                        collapsed: false
                        Layout.fillWidth: true
                        side_editor: App.materials.material_edit.back_editor
                    }
                }
            }
        }

        RowLayout {
            DeleteItemButton {
                id: delete_button

                title: qsTr("Delete Material")
                itemType: qsTr("material")
                replacementModel: root.source_model

                toDelete: AppData.materials.current_material

                onBeforeOpened: {
                    isDangerous = AppData.current_database.material_use_count(
                                delete_button.toDelete
                                )
                }

                onDeleteRequested: {
                    AppData.current_database.delete_material_group(
                                delete_button.toDelete
                                )
                    root.clearSelection()
                }

                onDeleteReassignRequested: (entity) => {
                    AppData.current_database.delete_material_group(
                                delete_button.toDelete,
                                entity
                                )
                    root.clearSelection()
                }
            }
        }
    }

    placeholder: Item {
        Label {
            anchors.centerIn: parent
            text: qsTr("Select a material")
            font.pointSize: 16
            opacity: 0.5
        }
    }
}
