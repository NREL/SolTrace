import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Material

import SolTrace

ListView {
    id: root
    width: 300
    height: 400
    clip: true

    ScrollBar.vertical: STScrollBar { }

    model: App.sun.shape.custom_distribution

    function nextAngle() {
        const count = root.model.count()
        if (count <= 0) {
            return 0
        }

        const lastAngle = root.model.data(root.model.index(count - 1, 0))
        return Math.min(20, lastAngle + 1)
    }

    header: RowLayout {
        width: root.width

        Label {
            text: qsTr("Angle (mrad)")
            Layout.fillWidth: true
            font.bold: true

            
        }

        Label {
            text: qsTr("Intensity")
            Layout.fillWidth: true
            font.bold: true
        }
    }

    delegate: RowLayout {
        width: root.width
        spacing: 4
        visible: index < root.model.rowCount()

        Label {
            id: index_label

            Layout.preferredWidth: 25
            text: index + 1

            
        }

        STDoubleSpinBox {
            id: angle

            from: 0
            to: 20
            decimals: 3
            Layout.fillWidth: true
            value: root.model.data(root.model.index(index, 0))
            onValueModified: root.model.setData(root.model.index(index, 0), value, Qt.EditRole)

            
        }

        STDoubleSpinBox {
            id: intensity

            from: 0
            to: 1.2
            stepSize: 0.1
            decimals: 3
            Layout.fillWidth: true
            value: root.model.data(root.model.index(index, 1))
            onValueModified: root.model.setData(root.model.index(index, 1), value, Qt.EditRole)
        }

        STIconButton {
            id: remove_button

            Layout.preferredWidth: implicitWidth
            Layout.preferredHeight: implicitWidth
            icon: "\uf00d"
            onClicked: if (App.sun.shape.custom_distribution.count() > 3) App.sun.shape.custom_distribution.remove(index)
            iconSize: 10
        }
    }

    footer: STButton {
        width: root.width

        text: qsTr("+ Add row")
        onClicked: App.sun.shape.custom_distribution.append(root.nextAngle(), 0)
    }
}
