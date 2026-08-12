import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property var module: App.layout.instance_edit
    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property var labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    STPropertyPanel {
        Layout.fillWidth: true
        columns: 2

        Label {
            Layout.columnSpan: 2
            Layout.fillWidth: true

            visible: root.module.current_material_name.length === 0
                     || root.module.current_geometry_name.length === 0
            text: "A material and geometry must be assigned for the element to be visible."
            color: Material.color(Material.Yellow)
            wrapMode: Label.WrapAtWordBoundaryOrAnywhere
        }

        STPropertyLabel {
            text: "Parent"
            Layout.alignment: root.labelAlignment
            Layout.columnSpan: root.singleColumn ? 2 : 1
        }

        STButton {
            Layout.fillWidth: true
            Layout.columnSpan: root.singleColumn ? 2 : 1
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

        STPropertyLabel {
            text: "Material"
            Layout.alignment: root.labelAlignment
            Layout.columnSpan: root.singleColumn ? 2 : 1
        }

        STButton {
            Layout.fillWidth: true
            Layout.columnSpan: root.singleColumn ? 2 : 1
            property string material_name: root.module.current_material_name
            text: material_name.length ? material_name : "Unassigned"
            onClicked: material_pop.open()
            SelectItemPopup {
                id: material_pop
                source_model: AppData.materials.materials_list
                onSelectedEntity: (entity) => root.module.current_material = entity
            }
        }

        STPropertyLabel {
            text: "Geometry"
            Layout.alignment: root.labelAlignment
            Layout.columnSpan: root.singleColumn ? 2 : 1
        }

        STButton {
            Layout.fillWidth: true
            Layout.columnSpan: root.singleColumn ? 2 : 1
            property string geometry_name: root.module.current_geometry_name
            text: geometry_name.length ? geometry_name : "Unassigned"
            onClicked: geometry_pop.open()
            SelectItemPopup {
                id: geometry_pop
                source_model: AppData.materials.geometry_list
                onSelectedEntity: (entity) => root.module.current_geometry = entity
            }
        }

        STPropertySeparator {
            title: "Placement"
        }

        STSwitch {
            Layout.columnSpan: 2
            text: "Use 3D widget"
            checked: App.view.mouse_mode === ViewModule.EditElement
            onToggled: App.view.mouse_mode = checked ? ViewModule.EditElement
                                                      : ViewModule.Camera
        }

        InlineDocumentation {
            key: "configure.layout.coordinates"
            Layout.columnSpan: 2
            Layout.fillWidth: true
        }

        STPropertyLabel {
            text: "Coordinates"
            Layout.alignment: root.labelAlignment
            Layout.columnSpan: root.singleColumn ? 2 : 1
        }

        STComboBox {
            id: positionModeCombo
            Layout.columnSpan: root.singleColumn ? 2 : 1
            Layout.fillWidth: true
            model: ["Local", "Global"]
        }

        StackLayout {
            id: positionSwipe
            Layout.columnSpan: 2
            Layout.fillWidth: true
            //Layout.preferredHeight: currentItem ? currentItem.implicitHeight : 0
            currentIndex: positionModeCombo.currentIndex

            GridLayout {
                id: localPositionEditor
                Layout.preferredWidth: positionSwipe.width
                columns: 2


                // Splitting out the positions to try and make this more robust
                property real xPosition: 0
                property real yPosition: 0
                property real zPosition: 0

                function set_position_values(x, y, z) {
                    xPosition = x
                    yPosition = y
                    zPosition = z
                    local_x_pos.value = x
                    local_y_pos.value = y
                    local_z_pos.value = z
                }

                function refresh_position() {
                    const position = root.module.position
                    set_position_values(position.x, position.y, position.z)
                }

                function update_position() {
                    root.module.position = Qt.vector3d(xPosition,
                                                       yPosition,
                                                       zPosition)
                }

                Component.onCompleted: refresh_position()

                Connections {
                    target: root.module

                    function onPosition_changed() {
                        localPositionEditor.refresh_position()
                    }
                }

                STPropertyLabel { text: "X" }
                STDoubleSpinBox {
                    id: local_x_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        localPositionEditor.xPosition = value
                        localPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STPropertyLabel { text: "Y" }
                STDoubleSpinBox {
                    id: local_y_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        localPositionEditor.yPosition = value
                        localPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STPropertyLabel { text: "Z" }
                STDoubleSpinBox {
                    id: local_z_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        localPositionEditor.zPosition = value
                        localPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STButton {
                    Layout.fillWidth: true
                    Layout.columnSpan: 2
                    text: "Reset Local Default Position"
                    onClicked: {
                        root.module.position = Qt.vector3d(0, 0, 0)
                        localPositionEditor.refresh_position()
                    }
                }
            }

            GridLayout {
                id: globalPositionEditor
                Layout.preferredWidth: positionSwipe.width
                columns: 2

                property real xPosition: 0
                property real yPosition: 0
                property real zPosition: 0

                function set_position_values(x, y, z) {
                    xPosition = x
                    yPosition = y
                    zPosition = z
                    global_x_pos.value = x
                    global_y_pos.value = y
                    global_z_pos.value = z
                }

                function refresh_position() {
                    const position = root.module.global_position
                    set_position_values(position.x, position.y, position.z)
                }

                function update_position() {
                    root.module.global_position =
                            Qt.vector3d(xPosition,
                                        yPosition,
                                        zPosition)
                }

                Component.onCompleted: refresh_position()

                Connections {
                    target: root.module

                    function onGlobal_position_changed() {
                        globalPositionEditor.refresh_position()
                    }
                }

                STPropertyLabel { text: "X" }
                STDoubleSpinBox {
                    id: global_x_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        globalPositionEditor.xPosition = value
                        globalPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STPropertyLabel { text: "Y" }
                STDoubleSpinBox {
                    id: global_y_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        globalPositionEditor.yPosition = value
                        globalPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STPropertyLabel { text: "Z" }
                STDoubleSpinBox {
                    id: global_z_pos
                    Layout.fillWidth: true
                    from: -Infinity
                    to: Infinity
                    onValueModified: {
                        globalPositionEditor.zPosition = value
                        globalPositionEditor.update_position()
                    }
                    decimals: 4
                }

                STButton {
                    Layout.fillWidth: true
                    Layout.columnSpan: 2
                    text: "Reset Global Default Position"
                    onClicked: {
                        root.module.global_position = Qt.vector3d(0, 0, 0)
                        globalPositionEditor.refresh_position()
                    }
                }
            }
        }

        STPropertySeparator {
            id: rotPanel
            title: "Parent-relative Rotation"

            property real xAngle: 0
            property real yAngle: 0
            property real zAngle: 0

            function set_angle_values(x, y, z) {
                xAngle = x
                yAngle = y
                zAngle = z
                x_euler.value = x
                y_euler.value = y
                z_euler.value = z
            }

            function refresh_angles() {
                const angles = root.module.euler_angles_xyz
                set_angle_values(angles.x, angles.y, angles.z)
            }

            function update_from_angles() {
                root.module.euler_angles_xyz = Qt.vector3d(xAngle,
                                                           yAngle,
                                                           zAngle)
            }

            Component.onCompleted: refresh_angles()

            Connections {
                target: root.module

                function onEuler_angles_xyz_changed() {
                    rotPanel.refresh_angles()
                }
            }
        }

            STPropertyLabel { text: "X Angle (deg)" }
            STDoubleSpinBox {
                id: x_euler
                Layout.fillWidth: true
                from: -Infinity
                to: Infinity
                onValueModified: {
                    rotPanel.xAngle = value
                    rotPanel.update_from_angles()
                }
            }

            STPropertyLabel { text: "Y Angle (deg)" }
            STDoubleSpinBox {
                id: y_euler
                Layout.fillWidth: true
                from: -Infinity
                to: Infinity
                onValueModified: {
                    rotPanel.yAngle = value
                    rotPanel.update_from_angles()
                }
            }

            STPropertyLabel { text: "Z Angle (deg)" }
            STDoubleSpinBox {
                id: z_euler
                Layout.fillWidth: true
                from: -Infinity
                to: Infinity
                onValueModified: {
                    rotPanel.zAngle = value
                    rotPanel.update_from_angles()
                }
            }

            STButton {
                Layout.fillWidth: true
                Layout.columnSpan: 2
                text: "Reset"
                onClicked: root.module.euler_angles_xyz = Qt.vector3d(0, 0, 0)
            }

            STButton {
                Layout.columnSpan: 2
                Layout.fillWidth: true
                text: "Point at..."
                onClicked: look_at_pop.open()

                STPopup {
                    id: look_at_pop

                    function accept_position() {
                        root.module.look_at_world_position(
                                    Qt.vector3d(look_at_x.value,
                                                look_at_y.value,
                                                look_at_z.value))
                    }

                    ColumnLayout {
                        spacing: 8

                        STComboBar {
                            id: look_at_mode
                            Layout.fillWidth: true
                            model: ["Position", "Element"]
                            iconModel: ["\uf3c5", "\uf6d1"]
                        }

                        GridLayout {
                            visible: look_at_mode.currentIndex === 0
                            Layout.fillWidth: true
                            columns: 2

                            STPropertyLabel { text: "X" }
                            STDoubleSpinBox {
                                id: look_at_x
                                Layout.fillWidth: true
                                value: 0
                                from: -Infinity
                                to: Infinity
                                decimals: 4
                                onValueModified: look_at_pop.accept_position()
                            }

                            STPropertyLabel { text: "Y" }
                            STDoubleSpinBox {
                                id: look_at_y
                                Layout.fillWidth: true
                                value: 0
                                from: -Infinity
                                to: Infinity
                                decimals: 4
                                onValueModified: look_at_pop.accept_position()
                            }

                            STPropertyLabel { text: "Z" }
                            STDoubleSpinBox {
                                id: look_at_z
                                Layout.fillWidth: true
                                value: 0
                                from: -Infinity
                                to: Infinity
                                decimals: 4
                                onValueModified: look_at_pop.accept_position()
                            }
                        }

                        STButton {
                            visible: look_at_mode.currentIndex === 0
                            Layout.fillWidth: true
                            text: "Point at Position"
                            onClicked: {
                                look_at_pop.accept_position()
                                look_at_pop.close()
                            }
                        }

                        STButton {
                            visible: look_at_mode.currentIndex === 1
                            Layout.fillWidth: true
                            text: "Choose Element"
                            onClicked: look_at_element_pop.open()

                            SelectElementPopup {
                                id: look_at_element_pop
                                exclude: [root.module.entity]
                                onSelectedElement: (element) => {
                                                      root.module.look_at_entity(element)
                                                      look_at_pop.close()
                                                  }
                            }
                        }
                    }
                }
            }

        STPropertySeparator {
            title: "Visualization"
        }

        STSwitch {
            Layout.columnSpan: 2
            text: "Hidden"
            checked: root.module.hidden
            onToggled: root.module.hidden = checked
        }

        ColorPickerField {
            id: elementColorPicker
            Layout.columnSpan: 2
            Layout.preferredWidth: 200
            color: root.module.color
            label: "Element Color"
            onUpdated: {
                root.module.color = elementColorPicker.color
            }
        }

        STSwitch {
            Layout.columnSpan: 2
            text: "Disabled"
            checked: root.module.disabled
            onToggled: root.module.disabled = checked
        }

        STSwitch {
            Layout.columnSpan: 2
            text: "Virtual"
            checked: root.module.virtual_element
            onToggled: root.module.virtual_element = checked
        }
    }
}
