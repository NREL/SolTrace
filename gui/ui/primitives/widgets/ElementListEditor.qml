import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

Item {
    id: root

    property alias model: internal.model

    property bool editing: false

    // Width threshold for switching between stacked and side-by-side
    property int wideThreshold: 500

    // List side preferred width in wide modes
    property int listWidth: 250

    // Header component shown above the list (e.g. a search field)
    property Component listHeader: null

    // Footer component shown below the list (e.g. an add button)
    property Component listFooter: null

    // The component shown when an element is selected
    property Component detailView: null

    // Placeholder shown in wide mode when nothing is selected
    property Component placeholder: null

    // Text shown over the list when the list model is empty
    property string emptyListText: ""

    // Read only data
    readonly property bool wideMode: width >= wideThreshold
    readonly property bool hasSelection: editing

    // Signal
    signal itemClicked(db_entity entity)

    onEditingChanged: {
        if (root.editing && narrowStack.depth < 2)
            narrowStack.push(narrowDetailComponent)
        else if (!root.editing && narrowStack.depth > 1)
            narrowStack.pop()
    }

    function goBack() {
        editing = false
    }

    function clearSelection() {
        editing = false
    }

    function _getModelCount() {
        if (!internal.model) {
            return 0
        }

        return (typeof internal.model.rowCount === "function")
                ? internal.model.rowCount()
                : internal.model.count
    }

    Connections {
        target: internal.model
        ignoreUnknownSignals: true
        function onRowsRemoved(modelParent, first, last) {
            root.editing = false
        }
    }

    QtObject {
        id: internal
        property var model: null
    }

    // Narrow stack layout
    Item {
        anchors.fill: parent
        visible: !root.wideMode

        StackView {
            id: narrowStack
            initialItem: narrowListComponent
            anchors.fill: parent
        }
    }

    Component{
        id: listDelegate

        STItemDelegate {
            id: st_del_root
            required property int index

            required property string name
            required property var entity

            signal opened(db_entity entity)

            text: name

            highlighted: entity === AppData.layout.edited_element

            onClicked: {
                st_del_root.ListView.view.editorRoot.editing = true
                st_del_root.ListView.view.editorRoot.itemClicked(entity)
            }

            onOpened: {
                st_del_root.ListView.view.editorRoot.editing = false
            }

            STIconButton {
                text: "\uf802"

                onClicked: {
                    AppData.layout.viewed_element = st_del_root.entity
                    AppData.layout.clear_edited_element()
                    st_del_root.opened(st_del_root.entity)
                }

                anchors.right: parent.right
                anchors.verticalCenter: parent.verticalCenter
            }
        }
    }

    Component {
        id: narrowListComponent

        ColumnLayout {
            spacing: 0

            Breadcrumb {
                Layout.fillWidth: true
                Layout.margins: 3

                onItemSelected: {
                    root.editing = true
                }

                onItemDeselected: {
                    root.editing = false
                }
            }

            Loader {
                Layout.fillWidth: true
                Layout.margins: 8

                sourceComponent: root.listHeader
                active: root.listHeader !== null
                visible: active
            }

            ListView {
                id: narrowListView
                Layout.fillHeight: true
                Layout.fillWidth: true
                clip: true
                model: internal.model
                property var editorRoot: root
                ScrollIndicator.vertical: ScrollIndicator { }

                delegate: listDelegate

                Label {
                    anchors.centerIn: parent
                    width: parent.width
                    horizontalAlignment: Text.AlignHCenter
                    wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                    text: root.emptyListText
                    opacity: 0.5
                    visible: root.emptyListText.length > 0
                             && narrowListView.count === 0
                }
            }

            Loader {
                Layout.fillWidth: true
                Layout.margins: 8
                sourceComponent: root.listFooter
                active: root.listFooter !== null
                visible: active
            }
        }
    }

    Component {
        id: narrowDetailComponent

        Loader {
            sourceComponent: root.detailView
            active: root.detailView !== null
        }
    }

   // Wide layout
    Item {
        anchors.fill: parent
        visible: root.wideMode

        RowLayout {
            anchors.fill: parent
            spacing: 0
            Layout.margins: 8


            ColumnLayout {
                Layout.preferredWidth: root.listWidth
                Layout.maximumWidth: root.listWidth
                Layout.fillHeight: true
                Layout.margins: 8
                spacing: 0

                Breadcrumb {
                    Layout.fillWidth: true
                    Layout.margins: 3

                    onItemSelected: {
                        root.editing = true
                    }

                    onItemDeselected: {
                        root.editing = false
                    }
                }

                Loader {
                    Layout.fillWidth: true
                    sourceComponent: root.listHeader
                    active: root.listHeader !== null
                    visible: active
                }

                ListView {
                    id: wideListView
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: internal.model
                    property var editorRoot: root
                    ScrollIndicator.vertical: ScrollIndicator { }

                    delegate: listDelegate

                    Label {
                        anchors.centerIn: parent
                        width: parent.width
                        horizontalAlignment: Text.AlignHCenter
                        wrapMode: Label.WrapAtWordBoundaryOrAnywhere
                        text: root.emptyListText
                        opacity: 0.5
                        visible: root.emptyListText.length > 0
                                 && wideListView.count === 0
                    }
                }

                Loader {
                    Layout.fillWidth: true
                    sourceComponent: root.listFooter
                    active: root.listFooter !== null
                    visible: active
                }
            }

            Rectangle {
                Layout.preferredWidth: 1
                Layout.fillHeight: true
                color: App.theme.dividerColor
            }

            Loader {
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.margins: 8
                sourceComponent: root.hasSelection ? root.detailView : root.placeholder
                active: true
            }
        }
    }
}
