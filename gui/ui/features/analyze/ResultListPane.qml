import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root

    property int selected_result_index: -1
    readonly property bool has_current_result: selected_result_index >= 0
                                               && selected_result_index < result_list.count

    signal closeRequested()

    spacing: 8

    function formatRayCount(count) {
        return Number(count).toLocaleString(Qt.locale(), "f", 0) + " rays"
    }

    function showResultMode() {
        App.view.simulation_content_view = true
        App.view.workflow_phase = ViewModule.Analyze
    }

    function showSceneMode() {
        App.view.simulation_content_view = false
    }

    Label {
        text: "Current Result"
        font.bold: true
    }

    RowLayout {
        Layout.fillWidth: true
        enabled: root.has_current_result
        spacing: 8

        Label {
            Layout.alignment: Qt.AlignVCenter
            font.family: "Font Awesome 7 Free"
            font.pointSize: 16
            text: "\uf303"
        }

        STTextField {
            Layout.fillWidth: true
            text: root.has_current_result ?
                      AppData.simulation.current_simulation_result_name : ""

            onTextEdited: {
                if (root.has_current_result) {
                    AppData.simulation.rename_result(root.selected_result_index,
                                                     text)
                }
            }
        }
    }

    RowLayout {
        Layout.fillWidth: true

        STIconButton {
            enabled: root.has_current_result
            icon: "\uf2ed"
            toolTip: "Delete Result"
            onClicked: AppData.simulation.delete_result(root.selected_result_index)
        }

        Item {
            Layout.fillWidth: true
        }

        STIconButton {
            enabled: root.has_current_result
            icon: "\uf24d"
            label: "Create Scene"
            toolTip: "Create Scene from Result"

            onClicked: {
                AppData.simulation.duplicate_current_result_for_edit()
                root.showSceneMode()
                root.closeRequested()
            }
        }

        STIconButton {
            enabled: root.has_current_result
            icon: "\uf019"
            label: "Export"
            toolTip: "Export Result"
            onClicked: {
                AppData.simulation.select_result(root.selected_result_index)
                AppData.exporter.export_current()
            }
        }
    }

    Rectangle {
        color: Material.dividerColor
        Layout.preferredHeight: 1
        Layout.fillWidth: true
        Layout.leftMargin: 3
        Layout.rightMargin: 3
    }

    InlineDocumentation {
        key: "analyze.results"
        Layout.fillWidth: true
    }

    Label {
        Layout.fillWidth: true
        text: "Simulation Results"
        font.bold: true
    }

    Label {
        Layout.fillWidth: true
        visible: result_list.count === 0
        text: "No simulation results yet."
        opacity: 0.7
        horizontalAlignment: Text.AlignHCenter
    }

    ListView {
        id: result_list
        property int previousCount: 0

        Layout.fillHeight: true
        Layout.fillWidth: true
        clip: true
        model: AppData.simulation.results

        ScrollBar.vertical: STScrollBar { }

        onCountChanged: {
            if (count > previousCount) {
                currentIndex = count - 1
                root.selected_result_index = currentIndex
            } else if (currentIndex >= count) {
                currentIndex = count - 1
                root.selected_result_index = currentIndex
            }
            previousCount = count
        }

        delegate: ItemDelegate {
            id: resultDelegate
            required property int index

            required property string name
            required property var when
            required property var ray_count

            width: result_list.width
            implicitHeight: 64
            highlighted: ListView.isCurrent ? ListView.isCurrent : false
            hoverEnabled: true

            background: Rectangle {
                opacity: resultDelegate.enabled ? 1 : 0.3
                color: resultDelegate.down
                       ? Material.rippleColor
                       : resultDelegate.highlighted
                         ? Material.highlightedRippleColor
                         : resultDelegate.hovered
                           ? Qt.alpha(Material.highlightedRippleColor, 0.35)
                           : "transparent"
                radius: 14
            }

            contentItem: RowLayout {
                spacing: 10

                ColumnLayout {
                    Layout.fillWidth: true
                    Layout.alignment: Qt.AlignVCenter
                    spacing: 2

                    Label {
                        Layout.fillWidth: true
                        text: resultDelegate.name
                        font.bold: true
                        font.pointSize: 18
                        elide: Text.ElideRight
                    }

                    Label {
                        Layout.fillWidth: true
                        text: Qt.formatDateTime(resultDelegate.when,
                                                "yyyy-MM-dd hh:mm:ss")
                        opacity: 0.7
                        elide: Text.ElideRight
                    }
                }

                Label {
                    Layout.alignment: Qt.AlignVCenter | Qt.AlignRight
                    text: root.formatRayCount(resultDelegate.ray_count)
                    opacity: 0.85
                    horizontalAlignment: Text.AlignRight
                }
            }

            onClicked: {
                result_list.currentIndex = index
                root.selected_result_index = index
                AppData.simulation.select_result(index)
                root.showResultMode()
            }
        }
    }
}
