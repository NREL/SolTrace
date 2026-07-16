import QtQuick
import QtQuick.Layouts
import SolTrace

RowLayout {
    id: root
    required property SplitPanelData target
    required property SplitPanelData otherTarget
    required property int available_width
    required property bool is_right_panel

    Layout.minimumWidth: implicitWidth
    spacing: 5

    STIconButton {
        id: inline_docs_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitWidth

        icon: "\uf059"
        toolTip: (root.target.inline_docs ? "Disable" : "Enable") + " Inline Docs"
        visible: root.target.size >= SplitPanelData.Normal
        onClicked: root.target.inline_docs = !root.target.inline_docs
    }

    STIconButton {
        id: popup_opts_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitWidth

        icon: "\uf2d2"
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

                    icon: "\uf422"
                    toolTip: "Smaller"
                    visible: root.target.size != SplitPanelData.Small
                    onClicked: {
                        if (root.target.size === SplitPanelData.Full)
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

                    icon: "\uf424"
                    toolTip: "Larger"
                    visible: root.target.size != SplitPanelData.Full
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

                    icon: "\uf065"
                    toolTip: "Full Size"
                    visible: root.target.size != SplitPanelData.Full
                    onClicked: {
                        root.target.width = root.target.sizes[SplitPanelData.Full]
                        App.view.fit_panels(root.available_width, root.is_right_panel, false)
                        window_opts_pop.close()
                    }
                }
            }
        }
    }

    STIconButton {
        id: close_button

        Layout.preferredWidth: implicitWidth
        Layout.preferredHeight: implicitWidth

        icon: "\uf00d"
        toolTip: "Close Panel"
        onClicked: {
            root.target.visible = false
            window_opts_pop.close()
        }
    }
}
