pragma Singleton

import QtQuick
import SolTrace

QtObject {
    id: root

    readonly property Theme theme: Theme {}

    readonly property var db: AppData.current_database
    readonly property var file_source: AppData.file_source
    readonly property var sun: AppData.sun
    readonly property var layout: AppData.layout
    readonly property var view: AppData.view
    readonly property var materials: AppData.materials
    readonly property var simulation: AppData.simulation
    readonly property var docs: AppData.docs


}
