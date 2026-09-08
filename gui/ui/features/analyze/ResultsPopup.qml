import QtQuick

import SolTrace

STPopup {
    id: root

    enable_shadow: true

    ResultListPane {
        anchors.fill: parent
        onCloseRequested: root.close()
    }
}
