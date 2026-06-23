import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Item {
    id: root
    implicitWidth: 200
    implicitHeight: 32

    property alias text: tf.text
    property alias placeholderText: tf.placeholderText

    property alias validator: tf.validator
    property alias acceptableInput: tf.acceptableInput

    signal accepted()
    signal textEdited(string text)

    property string leftIcon
    property string rightIcon

    WellRectangle {
            anchors.fill: parent
    }

    RowLayout {
        anchors.fill: parent

        spacing: 0

        Label {
            text: root.leftIcon

            Layout.fillHeight: true
            Layout.preferredWidth: height

            Layout.topMargin: 7
            Layout.bottomMargin: 7
            Layout.leftMargin: 8

            font.family: "Font Awesome 7 Free"

            font.pointSize: 48

            fontSizeMode: Label.Fit
            

            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            visible: text.length > 0
        }

        TextField {
            id: tf

            Layout.fillWidth: true
            Layout.fillHeight: true

            // These dont seem to work.
            leftInset: 0
            padding: 0

            background: Item {}
            

            onAccepted: root.accepted()
            onTextEdited: root.textEdited(text)

            // Rectangle {
            //     color: "red"
            //     anchors.fill: parent
            // }
        }

        Label {
            text: root.rightIcon

            Layout.fillHeight: true
            Layout.preferredWidth: height

            Layout.topMargin: 7
            Layout.bottomMargin: 7
            Layout.rightMargin: 8

            font.family: "Font Awesome 7 Free"

            font.pointSize: 48
            

            fontSizeMode: Label.Fit

            horizontalAlignment: Qt.AlignHCenter
            verticalAlignment: Qt.AlignVCenter

            visible: text.length > 0
        }
    }


}
