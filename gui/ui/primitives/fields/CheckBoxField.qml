import QtQuick
import SolTrace

STButton {
    property bool value: false

    text_icon: value ? "\uf205" : "\uf204"

    onClicked: {
        value = !value
    }
}
