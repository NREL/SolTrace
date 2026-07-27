import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts
import QtQuick.Dialogs
import SolTrace

ElementListEditor {
    id: root

    property var module: App.layout

    property bool has_viewed_entity: module.viewed_element.is_valid()

    model: has_viewed_entity ? module.filtered_child_model : module.filtered_root_elements_model
    wideThreshold: 500
    listWidth: 250
    emptyListText: "Scene elements combine geometry, material, and placement to define a physical object. Create one below to build the scene."

    onEditingChanged: {
        App.view.editing_layout = editing
        if (!editing && App.view.mouse_mode === ViewModule.EditElement) {
            App.view.mouse_mode = ViewModule.Camera
        }
    }

    Connections {
        target: App.view
        function onEditing_layout_changed() {
            root.editing = App.view.editing_layout
        }
    }

    onItemClicked: function(entity) {
        module.edited_element = entity
    }

    listHeader: RowLayout {
        STSearchField {
            id: search_field
            Layout.fillWidth: true
            text: root.model.name_filter

            Binding {
                target: root.model
                property: "name_filter"
                value: search_field.text
            }
        }

        STIconButton {
            icon: "\uf05b"
            toolTip: "Select Element From View"
            onClicked: App.view.mouse_mode = ViewModule.SelectElement
        }

        STIconButton {
            icon: "\uf0b0"

            onClicked: filter_popup.open()

            STPopup {
                id: filter_popup

                GridLayout {
                    columns: 3

                    Label {
                        Layout.fillWidth: true
                        Layout.columnSpan: 3

                        text: "Filter elements by:"
                    }

                    Label {
                        text: "Material:"
                    }

                    STClickableLabel {
                        Layout.fillWidth: true
                        property string mat_name: root.model.material_filter_name
                        text: mat_name.length ? mat_name : "All"

                        onClicked: mat_pop.open()

                        SelectItemPopup {
                            id: mat_pop
                            source_model: AppData.materials.materials_list

                            onSelectedEntity: (entity) =>
                                              root.model.material_filter = entity
                        }
                    }

                    STIconButton {
                        icon: "\uf0e2"
                        onClicked: root.model.clear_material()
                    }

                    Label {
                        text: "Geometry:"
                    }

                    STClickableLabel {
                        Layout.fillWidth: true
                        property string geo_name: root.model.geometry_filter_name
                        text: geo_name.length ? geo_name : "All"

                        onClicked: geo_pop.open()

                        SelectItemPopup {
                            id: geo_pop
                            source_model: AppData.materials.geometry_list

                            onSelectedEntity: (entity) =>
                                              root.model.geometry_filter = entity
                        }
                    }

                    STIconButton {
                        icon: "\uf0e2"
                        onClicked: root.model.clear_geometry()
                    }

                    STButton {
                        text: "Reset All"

                        onClicked: root.model.clear_all_filters()
                    }

                }
            }
        }

        STIconButton {
            id: clear_filter_button

            visible: root.model.has_filter

            icon: "\ue17b"

            onClicked: root.model.clear_all_filters()
        }
    }

    listFooter: RowLayout {
        CreateNewItemButton {
            title: "New Element"
            icon: "\uf055"

            onCreateRequested: function(name) {
                var new_name = AppData.current_database.sanitize_element_name(name)
                var entity = AppData.current_database.add_element(new_name)
                root.module.edited_element = entity
                root.editing = true
            }
        }

        CreateNewItemButton {
            title: "New Child Element"
            icon: "\uf0fe"
            label: "New Child"
            enabled: root.has_viewed_entity
            opacity: enabled ? 1.0 : 0.4

            onCreateRequested: function(name) {
                var new_name = AppData.current_database.sanitize_element_name(name)
                var entity = AppData.current_database.add_element(
                            new_name,
                            root.module.viewed_element
                            )
                root.module.edited_element = entity
                root.editing = true
            }
        }
    }

    detailView: ColumnLayout {

        RowLayout {
            STIconButton {
                icon: "\uf053"
                visible: !root.wideMode
                onClicked: {
                    root.module.clear_edited_element()
                    root.goBack()
                }
            }

            RenameLabel {
                text: root.module.edited_element_name

                onAccepted: (new_name) => {
                                var curr_db = root.module.current_database

                                new_name = curr_db.sanitize_element_name(new_name);

                                curr_db.set_name_of(
                                    root.module.edited_element,
                                    new_name
                                    )
                            }
            }
        }

        ScrollView {
            id: layout_scroll
            Layout.fillHeight: true
            Layout.fillWidth: true
            contentWidth: availableWidth

            ColumnLayout {
                width: layout_scroll.availableWidth

                InlineDocumentation {
                    key: "configure.layout"
                }

                ElementEdit {
                    Layout.fillWidth: true
                }
            }
        }

        RowLayout {


            STIconButton {
                icon: "\uf2ed"
                label: "Delete"
                onClicked: delete_element_dialog.open()

                STDialog {
                    id: delete_element_dialog

                    title: "Delete Element"

                    ColumnLayout {
                        anchors.fill: parent

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.WordWrap
                            text: "Delete this element? Any children of this element will be made top-level."
                        }
                    }

                    footer: STDialogButtonBox {
                        STButton {
                            text: qsTr("Cancel")
                            DialogButtonBox.buttonRole: DialogButtonBox.RejectRole
                        }

                        STButton {
                            text: qsTr("Delete")
                            Material.foreground: Material.Red
                            DialogButtonBox.buttonRole: DialogButtonBox.DestructiveRole
                            idle_color: App.theme.destructiveGlassColor
                            down_color: Material.color(Material.Red)

                            onClicked: {
                                root.module.delete_edited_element()
                                root.goBack()
                                delete_element_dialog.close()
                            }
                        }
                    }
                }
            }

            Item {
                Layout.fillWidth: true
            }

            STIconButton {
                icon: "\uf140"
                label: "Orient Camera"
                toolTip: "Orient Camera to Element"
                enabled: root.module.instance_edit
                onClicked: simulation_scene.orient_camera_to_database_position(
                               root.module.instance_edit.global_position)
            }
        }
    }

    placeholder: Item {
        Label {
            anchors.centerIn: parent
            text: "Select an element"
            font.pointSize: 16
            opacity: 0.5
        }
    }
}
