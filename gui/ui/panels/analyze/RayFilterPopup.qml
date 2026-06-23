import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

STPopup {
    id: root

    signal modified(list<string> filter)

    property bool updating: false

    margins: 1
    padding: 8

    function open_with(initial_filter) {
        const filter = initial_filter || []
        updating = true

        for (let i = 0; i < ray_type_repeater.count; ++i) {
            const item = ray_type_repeater.itemAt(i)
            if (item) {
                item.checked = filter.indexOf(item.ray_type) !== -1
            }
        }

        updating = false

        open()
    }

    function current_filter() {
        let filter = []

        for (let i = 0; i < ray_type_repeater.count; ++i) {
            const item = ray_type_repeater.itemAt(i)
            if (item && item.checked) {
                filter.push(item.ray_type)
            }
        }

        return filter
    }

    property var ray_type_model: [
        { label: "Create", ray_type: "create" },
        { label: "Absorb", ray_type: "absorb" },
        { label: "Reflect", ray_type: "reflect" },
        { label: "Transmit", ray_type: "transmit" },
        { label: "Virtual", ray_type: "virtual" },
        { label: "Exit", ray_type: "exit" },
    ]

    GridLayout {
        columns: 2

        ColumnLayout {
            Repeater {
                id: ray_type_repeater
                model: ray_type_model
                delegate: STSwitch {
                    Layout.fillWidth: true
                    property string ray_type: modelData.ray_type
                    text: modelData.label
                    onToggled: {
                        if (!root.updating) {
                            root.modified(root.current_filter())
                        }
                    }
                }
            }
        }
    }
}
