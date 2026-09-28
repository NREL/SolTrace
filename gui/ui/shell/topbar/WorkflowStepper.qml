import QtQuick
import QtQuick.Layouts
import SolTrace

RowLayout {
    id: root

    property int currentIndex: 0
    property string previous: ""
    property string next: ""

    STButton {
        left_text_icon: "\uf177"
        text: root.previous
        visible: root.currentIndex > ViewModule.Load
        onClicked: {
            App.view.workflow_phase = root.currentIndex - 1
        }
    }

    Item {
        Layout.fillWidth: true
    }

    STButton {
        right_text_icon: "\uf178"
        text: root.next
        visible: root.currentIndex < ViewModule.Analyze
        onClicked: {
            App.view.workflow_phase = root.currentIndex + 1
        }
    }
}
