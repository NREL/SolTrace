import QtQuick
import QtQuick.Layouts
import SolTrace

RowLayout {
    id: root
    required property PanelData target
    required property PanelData otherTarget
    required property int available_width
    required property bool is_right_panel

    Layout.minimumWidth: implicitWidth
    spacing: 5

    STIconButton {
        id: inline_docs_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitWidth

        text: "\uf05a"
        toolTip: (root.target.inline_docs ? "Disable" : "Enable") + " Inline Docs"
        visible: root.target.size >= PanelData.Normal
        onClicked: root.target.inline_docs = !root.target.inline_docs
    }

    STIconButton {
        id: popup_opts_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitWidth

        text: "\uf2d2"
        toolTip: "Resize Panel"
        onClicked: window_opts_pop.open()

        STPopup {
            id: window_opts_pop
            RowLayout {
                anchors.fill: parent

                STIconButton {
                    id: smaller_button

                    Layout.preferredWidth: implicitWidth
                    Layout.preferredHeight: implicitWidth

                    text: "\uf422"
                    toolTip: "Smaller"
                    visible: root.target.size != PanelData.Small
                    onClicked: {
                        if (root.target.size === PanelData.Full)
                            root.otherTarget.visible = true

                        if (root.target.size >= 1){
                            root.target.width = root.target.sizes[root.target.size - 1]
                        }
                        App.view.fit_panels(root.available_width, root.is_right_panel, false)
                        window_opts_pop.close()
                    }
                }

                STIconButton {
                    id: larger_button

                    Layout.preferredWidth: implicitWidth
                    Layout.preferredHeight: implicitWidth

                    text: "\uf424"
                    toolTip: "Larger"
                    visible: root.target.size != PanelData.Full
                    onClicked: {
                        root.target.width = root.target.sizes[root.target.size + 1]
                        App.view.fit_panels(root.available_width, root.is_right_panel, false)
                        window_opts_pop.close()
                    }
                }

                STIconButton {
                    id: fullsize

                    Layout.preferredWidth: implicitWidth
                    Layout.preferredHeight: implicitWidth

                    text: "\uf065"
                    toolTip: "Full Size"
                    visible: root.target.size != PanelData.Full
                    onClicked: {
                        root.target.width = root.target.sizes[PanelData.Full]
                        App.view.fit_panels(root.available_width, root.is_right_panel, false)
                        window_opts_pop.close()
                    }
                }

                STIconButton {
                    id: close_button

                    Layout.preferredWidth: implicitWidth
                    Layout.preferredHeight: implicitWidth

                    text: "\uf00d"
                    toolTip: "Close Panel"
                    onClicked: {
                        root.target.visible = false
                        window_opts_pop.close()
                    }
                }
            }
        }
    }
}
