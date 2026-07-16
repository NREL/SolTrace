import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material
import SolTrace

ColumnLayout {
    spacing: 8

    Label {
        Layout.fillWidth: true
        Layout.leftMargin: 10
        Layout.rightMargin: 10
        Layout.bottomMargin: 8

        text: "Add a scene to the environment by loading a .stinput file"
        wrapMode: Text.WordWrap
    }

    SceneListPane {
        Layout.fillWidth: true
        Layout.fillHeight: true
    }

    WorkflowStepper {
        previous: "Get Started"
        next: "Configure Scene"
        currentIndex: ViewModule.Load
    }
}
