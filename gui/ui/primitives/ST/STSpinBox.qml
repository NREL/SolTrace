import QtQuick
import QtQuick.Controls.Material
import SolTrace

SpinBox {
    id: control
    editable: true
    live: true

    property string suffix

    function commitText(restoreInvalid) {
        const nextValue = Math.max(
                    control.from,
                    Math.min(control.to,
                             control.valueFromText(input.text, control.locale)))

        if (Number.isNaN(nextValue)) {
            if (restoreInvalid) {
                input.text = control.textFromValue(control.value, control.locale)
            }
            return
        }

        if (nextValue === control.value) {
            input.text = control.textFromValue(control.value, control.locale)
            return
        }

        control.value = nextValue
        control.valueModified()
    }

    contentItem: TextInput {
        id: input
        readonly property int suffixRightMargin: 30
        readonly property int suffixSpacing: 8

        z: 2
        color: App.theme.fontColor
        text: control.textFromValue(control.value, control.locale)
        font.family: control.font.family
        font.pointSize: App.theme.labelSize
        horizontalAlignment: Qt.AlignHCenter
        verticalAlignment: Qt.AlignVCenter
        rightPadding: suffixLabel.visible ? suffixLabel.width + suffixRightMargin + suffixSpacing : 0
        readOnly: !control.editable
        validator: control.validator
        inputMethodHints: control.inputMethodHints
        onTextEdited: {
            if (control.live && acceptableInput) {
                control.commitText(false)
            }
        }
        onAccepted: control.commitText(true)
        onEditingFinished: control.commitText(true)

        Label {
            id: suffixLabel
            anchors.right: parent.right
            anchors.rightMargin: input.suffixRightMargin
            anchors.verticalCenter: parent.verticalCenter
            color: App.theme.fontColor
            font.family: input.font.family
            font.pointSize: input.font.pointSize
            text: control.suffix
            visible: text.length > 0
        }
    }

    background: WellRectangle {
        implicitWidth: 80
        implicitHeight: 32
        radius: height / 2
    }

    up.indicator: Rectangle {
        x: control.mirrored ? 0 : parent.width - width
        height: parent.height
        implicitWidth: 32
        implicitHeight: 32
        color: control.up.pressed ? Qt.rgba(1, 1, 1, 0.1) : "transparent"
        radius: height / 2
        opacity: (control.hovered || control.activeFocus) ? 1 : 0

        Label {
            text: "\u002b"
            font.family: "Font Awesome 7 Free"
            font.pointSize: control.font.pointSize

            anchors.centerIn: parent
        }
    }

    down.indicator: Rectangle {
        x: control.mirrored ? parent.width - width : 0
        height: parent.height
        implicitWidth: 32
        implicitHeight: 32
        color: control.down.pressed ? Qt.rgba(1, 1, 1, 0.1) : "transparent"
        radius: height / 2
        opacity: (control.hovered || control.activeFocus) ? 1 : 0

        Label {
            text: "\uf068"
            font.family: "Font Awesome 7 Free"
            font.pointSize: control.font.pointSize

            anchors.centerIn: parent
        }
    }
}
