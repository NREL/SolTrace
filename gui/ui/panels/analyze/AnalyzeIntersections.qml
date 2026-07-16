import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

Flickable {
    id: root
    property var left_panel_size: App.view.left_panel.size
    property var intersections_module : AppData.intersections
    property bool singleColumn: App.view.left_panel.size === SplitPanelData.Small
    property int labelAlignment: (singleColumn ? Qt.AlignLeft : Qt.AlignRight) | Qt.AlignVCenter

    property var ray_geom: AppData.intersections.ray_geometry

    property int child_column_span: singleColumn ? 1 : 2

    property bool is_selected_ray_valid: ray_geom.selected_ray_id >= 0

    function formatRayCount(count) {
        return Number(count).toLocaleString(Qt.locale(), "f", 0) + " rays"
    }

    function visibleRayCount() {
        const available = ray_geom.available_rays
        return Math.round(available *
                          ray_geom.show_percent / 100.0)
    }

    function setVisibleRayCount(count) {
        const available = ray_geom.available_rays
        ray_geom.show_percent =
                available > 0 ? count * 100.0 / available : 0.0
    }

    contentWidth: width
    contentHeight: content_column.implicitHeight
    clip: true
    boundsBehavior: Flickable.StopAtBounds



    ColumnLayout {
        id: content_column
        width: root.width

        InlineDocumentation {
            key: "analyze.intersections"
            target: App.view.left_panel
        }

        Label {
            Layout.fillWidth: true
            Layout.columnSpan: 2
            text: root.formatRayCount(
                      root.ray_geom.available_rays)
            opacity: 0.75
            horizontalAlignment: Text.AlignHCenter
        }

        STPropertyPanel {
            Layout.fillWidth: true
            columns: root.singleColumn ? 1 : 2

            collapsible: true
            title: "Ray Visibility"


            STSwitch {
                Layout.fillWidth: true
                Layout.columnSpan: root.child_column_span

                text: "Show Intersections"
                checked: AppData.view.show_intersections

                onToggled: AppData.view.show_intersections = checked
            }

            STPropertyLabel {
                text: "Filter types"
                Layout.alignment: root.labelAlignment
            }

            STButton {
                Layout.fillWidth: true
                text: root.ray_geom.event_include.length > 0
                      ? root.ray_geom.event_include.join(", ")
                      : "None"

                onClicked: ray_filter_popup.open_with(
                               root.ray_geom.event_include)

                RayFilterPopup {
                    id: ray_filter_popup

                    onModified: function(filter) {
                        root.ray_geom.event_include = filter
                    }
                }
            }

            STDoubleSpinBox {
                Layout.columnSpan: root.child_column_span
                from: 0.0
                to: 100.0
                value: root.ray_geom.show_percent
                stepSize: 1.0
                Layout.fillWidth: true
                onValueModified: {
                    root.ray_geom.show_percent = value
                }
                suffix: "%"
            }

            STDoubleSpinBox {
                Layout.columnSpan: root.child_column_span
                from: 0.0
                to: root.ray_geom.available_rays
                value: root.visibleRayCount()
                stepSize: 1000
                decimals: 0
                Layout.fillWidth: true
                onValueModified: root.setVisibleRayCount(value)
                suffix: "rays"
            }

            STPropertySeparator {
                title: "Selection"
            }

            STPropertyLabel {
                text: "Selected ID"
                Layout.alignment: root.labelAlignment

                visible: root.is_selected_ray_valid
            }

            RowLayout {
                visible: root.is_selected_ray_valid

                Label {
                    Layout.fillWidth: true
                    text: root.ray_geom.selected_ray_id

                    visible: root.is_selected_ray_valid
                }

                STIconButton {
                    icon: "\uf057"

                    onClicked: root.ray_geom.selected_ray_id = -1
                }
            }

            STButton {
                Layout.columnSpan: root.child_column_span
                Layout.fillWidth: true
                text: "Pick Ray From View"
                left_text_icon: "\uf05b"
                onClicked: App.view.mouse_mode = ViewModule.PickRay
            }


        }
    }
}
