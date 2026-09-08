import QtCore
import QtQuick
import QtQuick.Layouts
import QtQuick.Dialogs
import QtQuick.Controls.Material

import SolTrace

ListView {
    id: root

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

        Item {
            Layout.preferredWidth: 25
        }

        Label {
            text: qsTr("Angle (mrad)")
            //Layout.fillWidth: true
            font.bold: true

            Layout.preferredWidth: parent.width * .4
        }

        Label {
            text: qsTr("Intensity")
            Layout.fillWidth: true
            font.bold: true
        }
    }


    delegate: RowLayout {
        id: row_del
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
            Layout.preferredWidth: row_del.width * .4
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

    // Not optimal, but for now
    footer: ColumnLayout{
        STButton {
            width: root.width

            text: qsTr("+ Add row")
            onClicked: App.sun.shape.custom_distribution.append(root.nextAngle(), 0)
        }
        RowLayout {
            STIconButton {
                icon: "\uf56f"
                label: "From CSV"

                onClicked: import_dialog.open()

                FileDialog {
                    id: import_dialog
                    nameFilters: ["CSV files (*.csv)"]
                    currentFolder: StandardPaths.standardLocations(StandardPaths.DocumentsLocation)[0]
                    onAccepted: App.sun.shape.import_custom_distribution(currentFile)
                }
            }

            STIconButton {
                icon: "\uf56e"
                label: "To CSV"

                onClicked: export_dialog.open()

                FileDialog {
                    id: export_dialog
                    nameFilters: ["CSV files (*.csv)"]
                    fileMode: FileDialog.SaveFile
                    currentFolder: StandardPaths.standardLocations(StandardPaths.DocumentsLocation)[0]
                    onAccepted: App.sun.shape.export_custom_distribution(currentFile)
                }
            }
        }
    }
}
