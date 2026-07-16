import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ColumnLayout {
    id: root
    required property string key

    property SplitPanelData target
    property string titlePrefix: ""
    property string title: ""
    property string body: ""

    visible: root.target ? root.target.inline_docs && root.target.size > 0 : true

    SubHeader {
        id: title
        text: {
            App.docs.locale
            App.docs.version
            return root.title ? root.title : App.docs.get(root.key, "title")
        }

        Layout.fillWidth: true
    }

    Label {
        id: body

        Layout.fillWidth: true

        textFormat: Label.MarkdownText

        text: {
            App.docs.locale
            App.docs.version
            return root.body ? root.body : App.docs.get(root.key)
        }

        wrapMode: Text.WordWrap

        font.pointSize: App.theme.normalSize
        color: App.theme.fontColor
    }
}
