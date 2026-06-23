import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

AdaptiveEditor {
    id: root

    property var source_model: null

    property string filterText

    onFilterTextChanged: {
        filtered_model.invalidate()
        clearSelection()
    }

    function selectSourceIndex(sourceIndex) {
        if (sourceIndex < 0) {
            clearSelection()
            return
        }

        try {
            var sourceIndexObject = root.source_model.index(sourceIndex, 0)
            if (filtered_model.mapFromSource) {
                var proxyIndexObject = filtered_model.mapFromSource(sourceIndexObject)
                if (proxyIndexObject && proxyIndexObject.row >= 0) {
                    currentIndex = proxyIndexObject.row
                    editing = true
                    return
                }

                clearSelection()
                return
            }
        } catch (error) {
            console.log("[AdaptiveFilteredEditor] source index mapping failed:",
                        error)
        }

        currentIndex = sourceIndex
        editing = true
    }

    model: SortFilterProxyModel {
        id: filtered_model
        model: root.source_model

        filters: [
            FunctionFilter {
                id: filter

                component NameFilter: QtObject {
                    property string name
                }

                function filter(data: NameFilter) : bool {
                    if (root.filterText.length) {
                        return (data.name || "").toLowerCase().includes(
                                    root.filterText.toLowerCase()
                                    )
                    }
                    return true
                }
            }
        ]
    }
}
