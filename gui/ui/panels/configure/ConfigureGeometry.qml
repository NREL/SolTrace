import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

AdaptiveFilteredEditor {
    id: root

    source_model: App.materials.geometry_list
    wideThreshold: 500
    listWidth: 250
    emptyListText: "Geometries define reusable element shapes. Create one below to assign it to elements."

    onEditingChanged: {
        App.view.editing_geometry = editing
    }

    Connections {
        target: App.view
        function onEditing_geometry_changed() {
            root.editing = App.view.editing_geometry
        }
    }

    Connections {
        target: App.materials
        function onCurrent_geometry_changed() {
            root.selectSourceIndex(
                        root.source_model.index_of(App.materials.current_geometry))
        }
    }

    onItemClicked: function(index, modelData) {
        App.materials.current_geometry = modelData.entity
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
            toolTip: "Select Geometry From View"
            onClicked: App.view.mouse_mode = ViewModule.SelectGeometry
        }
    }

    listFooter: RowLayout {
        CreateNewItemButton {
            title: "New Geometry"

            onCreateRequested: function(name) {
                var new_name = AppData.current_database.sanitize_geometry_name(name)
                AppData.current_database.add_geometry_group(new_name)
            }
        }
    }

    listDelegate: STItemDelegate {
        text: itemModel ? itemModel.name : "No name"
        highlighted: isCurrent
    }

    detailView: ColumnLayout {

        RowLayout {
            STIconButton {
                icon: "\uf053"
                visible: !root.wideMode
                onClicked: root.goBack()
            }

            RenameLabel {
                text: App.materials.current_geometry_name

                onAccepted: (new_name) => {
                    var mats = App.materials;
                    var curr_db = mats.current_database

                    new_name = curr_db.sanitize_geometry_name(new_name);

                    curr_db.set_name_of(
                                    mats.current_geometry,
                                    new_name
                                    )
                }
            }
        }

        SurfacePreviewScene {
            id: preview_scene
            Layout.fillWidth: true

            property bool show_expand: false

            Layout.preferredHeight: show_expand ? 300 : 148

            property bool no_geometry: App.materials.geometry_edit.surface_geometry.vertex_count === 0


            STIconButton {
                anchors.top: parent.top
                anchors.left: parent.left

                icon: "\uf424"

                onClicked: preview_scene.show_expand = !preview_scene.show_expand
            }

            // STIconButton {
            //     anchors.top: parent.top
            //     anchors.right: parent.right

            //     text: "\ue80b"

            //     onClicked: {

            //     }
            // }

            RowLayout {
                visible: parent.no_geometry

                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                anchors.margins: 3

                Label {
                    text: "\uf071"

                    font.family: "Font Awesome 7 Free"

                    color: Material.color(Material.Yellow)
                }

                Label {
                    Layout.fillWidth: true
                    text: "No geometry for this configuration."

                    elide: Label.ElideRight

                    color: Material.color(Material.Yellow)
                }
            }
        }

        ScrollView {
            id: geometry_scroll
            Layout.fillHeight: true
            Layout.fillWidth: true
            contentWidth: availableWidth

            GeometryProperties {
                width: geometry_scroll.availableWidth
                height: implicitHeight
            }
        }

        RowLayout {
            DeleteItemButton {
                id: delete_button

                title: "Delete Geometry"
                itemType: "geometry"
                replacementModel: root.source_model

                toDelete: AppData.materials.current_geometry

                onBeforeOpened: {
                    isDangerous = AppData.current_database.geometry_use_count(
                                delete_button.toDelete
                                )
                }

                onDeleteRequested: {
                    AppData.current_database.delete_geometry_group(
                                delete_button.toDelete
                                )
                    root.clearSelection()
                }

                onDeleteReassignRequested: (entity) => {
                    AppData.current_database.delete_geometry_group(
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
            text: "Select a geometry"
            font.pointSize: 16
            opacity: 0.5
            
        }
    }
}
