import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import SolTrace

STPropertyPanel {
    id: root
    property bool has_selection : false
    property var material_editor : App.materials.group_edit
    property var side_editor
    property bool singleColumn: App.view.left_panel.size <= SplitPanelData.Wide
    property var labelAlignment: (root.singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    columns: root.singleColumn ? 1 : 2

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Preset")
        Layout.alignment: root.labelAlignment
    }

    STButton {
        text: qsTr("Load")

        Layout.fillWidth: true

        onClicked: {
            preset_popup.open()
        }

        STPopup {
            id: preset_popup

            RowLayout {
                Layout.fillWidth: true

                STButton {
                    Layout.fillWidth: true
                    text: qsTr("Absorber")
                    onClicked: {
                        if (root.side_editor) {
                            root.side_editor.set_ideal_absorption()
                        }
                        preset_popup.close()
                    }
                }

                STButton {
                    Layout.fillWidth: true
                    text: qsTr("Reflector")
                    onClicked: {
                        if (root.side_editor) {
                            root.side_editor.set_ideal_reflection()
                        }
                        preset_popup.close()
                    }
                }

                STButton {
                    Layout.fillWidth: true
                    text: qsTr("Transmitter")
                    onClicked: {
                        if (root.side_editor) {
                            root.side_editor.set_ideal_transmission()
                        }
                        preset_popup.close()
                    }
                }
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Interaction")
        Layout.alignment: root.labelAlignment
    }

    STComboBox {
        id: interactionCombo
        Layout.fillWidth: true
        model: App.materials.material_edit.interaction_type_model
        textRole: "display"
        valueRole: "display"

        currentValue: root.side_editor.interaction_type

        displayText: currentText.charAt(0).toUpperCase() + currentText.slice(1).toLowerCase();

        onActivated: {
            if (root.side_editor && currentText.length > 0) {
                root.side_editor.interaction_type = currentText
            }
        }
    }



    // =========================================================================

    STPropertyLabel {
        text: qsTr("Distribution")
        Layout.alignment: root.labelAlignment
    }

    STComboBox {
        id: distributionCombo
        Layout.fillWidth: true
        model: App.materials.material_edit.distribution_type_model
        textRole: "display"
        valueRole: "display"

        currentValue: root.side_editor.error_distribution_type

        displayText: currentText.charAt(0).toUpperCase() + currentText.slice(1).toLowerCase();

        onActivated: {
            if (root.side_editor && currentText.length > 0) {
                root.side_editor.error_distribution_type = currentText
            }
        }
    }


    // =========================================================================

    STPropertyLabel {
        text: qsTr("Reflectance (ρ)")
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 0.0
        to: 1.0
        decimals: 4
        stepSize: 0.01
        value: root.side_editor ? root.side_editor.reflectivity : 0
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.reflectivity = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Transmittance (τ)")
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 0
        to: 1.0
        decimals: 4
        stepSize: 0.01
        value: root.side_editor ? root.side_editor.transmissivity : 0
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.transmissivity = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("<em>n</em> Front")
        textFormat: Label.RichText
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 1.0
        to: 5.0
        decimals: 4
        stepSize: 0.01
        value: root.side_editor ? root.side_editor.refraction_index_front : 1
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.refraction_index_front = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("<em>n</em> Back")
        textFormat: Label.RichText
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 1
        to: 5
        decimals: 4
        stepSize: 0.01
        value: root.side_editor ? root.side_editor.refraction_index_back : 1
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.refraction_index_back = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Slope Error (σ<sub>slope</sub>)")
        textFormat: Label.RichText
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 0
        to: 1000
        decimals: 3
        stepSize: 0.01
        suffix: "mrad"
        value: root.side_editor ? root.side_editor.slope_error : 0
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.slope_error = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Specularity Error (σ<sub>spec</sub>)")
        textFormat: Label.RichText
        Layout.alignment: root.labelAlignment
    }

    STDoubleSpinBox {
        Layout.fillWidth: true
        from: 0
        to: 1000
        decimals: 3
        stepSize: 0.01
        suffix: "mrad"
        value: root.side_editor ? root.side_editor.specularity_error : 0
        onValueModified: {
            if (root.side_editor) {
                root.side_editor.specularity_error = value
            }
        }
    }

    // =========================================================================

    STPropertyLabel {
        text: qsTr("Error Type")
        visible: false
        Layout.alignment: root.labelAlignment
    }

    STComboBox {
        Layout.fillWidth: true
        model: App.materials.material_edit.distribution_type_model
        textRole: "display"
        enabled: count > 0
        visible: false
        onActivated: {
            if (root.side_editor && currentText.length > 0) {
                root.side_editor.error_distribution_type = currentText
            }
        }
    }

    // =========================================================================

    // RowLayout {
    //     Layout.columnSpan: 2
    //     Layout.fillWidth: true
    //     spacing: 12

    //     Button {
    //         Layout.fillWidth: true
    //         text: "Angular Reflectance"
    //         enabled: false
    //     }

    //     Button {
    //         Layout.fillWidth: true
    //         text: "Angular Transmittance"
    //         enabled: false
    //     }
    // }
}
