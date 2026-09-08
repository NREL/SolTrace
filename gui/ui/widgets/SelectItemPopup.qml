import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

STPopup {
    id: root

    property var source_model: null

    property bool allowNothing: false

    property string filterText

    property list<db_entity> exclude

    signal selectedEntity(db_entity item)
    signal selectedNothing()

    onFilterTextChanged: {
        filtered_model.invalidate()
    }

    onExcludeChanged: {
        filtered_model.invalidate()
    }

    ColumnLayout {
        anchors.fill: parent

        spacing: 0

        STSearchField {
            id: search_field
            Layout.fillWidth: true

            Binding {
                target: root
                property: "filterText"
                value: search_field.text
            }
        }

        ListView {
            Layout.preferredHeight: 256
            Layout.fillWidth: true
            clip: true

            component NameEntityPair: QtObject {
                property string name
                property db_entity entity
            }

            model: SortFilterProxyModel {
                id: filtered_model
                model: root.source_model

                filters: [
                    FunctionFilter {
                        id: filter

                        function filter(data: NameEntityPair) : bool {
                            if (root.filterText.length) {
                                return (data.name || "").toLowerCase().includes(
                                            root.filterText.toLowerCase()
                                            )
                            }
                            return true
                        }
                    },
                    FunctionFilter {
                        id: exclude_filter

                        enabled: root.exclude.length > 0

                        inverted: true

                        function filter(data: NameEntityPair) : bool {
                            return root.exclude.includes(data.entity)
                        }
                    }
                ]
            }

            ScrollBar.vertical: STScrollBar { }

            delegate: STItemDelegate {
                id: delegate

                required property int index
                required property string name
                required property var entity

                text: name

                onClicked: {
                    selectedEntity(delegate.entity)
                    root.close()
                }
            }
        }

        STButton {
            id: do_nothing
            visible: root.allowNothing
            Layout.fillWidth: true

            text: qsTr("Unassign")

            onClicked: {
                root.selectedNothing()
                root.close()
            }
        }
    }
}
