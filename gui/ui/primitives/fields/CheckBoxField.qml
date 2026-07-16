import QtQuick
import SolTrace

STButton {
    property bool value: false

    left_text_icon: value ? "\uf205" : "\uf204"

    onClicked: {
        value = !value
    }
}
