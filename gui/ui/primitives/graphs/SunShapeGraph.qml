import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Material
import QtGraphs
import SolTrace

Rectangle {
    id: root

    Connections {
        target: App.sun.shape
        function onChanged() {
            rebuild()
        }
    }

    // Graph content
    property string title: ""
    property string xAxisTitle: ""
    property string yAxisTitle: ""

    // Axis ranges
    property real xMin: 0
    property real xMax: 1
    property real yMin: 0
    property real yMax: 1

    // Series access
    property alias series: lineSeries

    // Model
    property var model: App.sun.shape.current_distribution

    function rebuild() {
        lineSeries.clear()
        if (!model) {
            return;
        }

        const points = []
        for (let i = 0; i < model.rowCount(); i++) {
            const angle = Math.abs(model.data(model.index(i, 0)))
            const intensity = model.data(model.index(i, 1))
            points.push({
                            "angle": angle,
                            "intensity": intensity
                        })
        }

        points.sort(function(a, b) {
            return a.angle - b.angle
        })

        const mergedPoints = []
        for (let j = 0; j < points.length; ++j) {
            const previous = mergedPoints[mergedPoints.length - 1]
            if (previous && Math.abs(previous.angle - points[j].angle) < 1.0e-6) {
                previous.intensity = Math.max(previous.intensity, points[j].intensity)
            } else {
                mergedPoints.push(points[j])
            }
        }

        for (let k = mergedPoints.length - 1; k >= 0; --k) {
            if (mergedPoints[k].angle <= 1.0e-9) {
                continue
            }
            lineSeries.append(-mergedPoints[k].angle, mergedPoints[k].intensity)
        }

        for (let m = 0; m < mergedPoints.length; ++m) {
            lineSeries.append(mergedPoints[m].angle, mergedPoints[m].intensity)
        }
    }

    onModelChanged: {
        rebuild()
    }

    color: "transparent"
    radius: 8

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 10
        spacing: 0

        Item {
            id: graphHeader
            Layout.fillWidth: true
            Layout.preferredHeight: Math.max(graphTitle.implicitHeight, graphButtons.implicitHeight)

            Label {
                id: graphTitle
                anchors.centerIn: parent
                text: root.title
                font.bold: true
            }

            STIconButton {
                id: graphButtons
                anchors.right: parent.right
                text: "\uf1de"
                onClicked: graphSettings.open()

                STPopup {
                    id: graphSettings

                    ColumnLayout {
                        width: parent.width

                        STSpinBoxField {
                            label: "Line Width"
                            Layout.fillWidth: true
                            value: App.theme.sunShapeGraphLineWidth
                            onValueChanged: App.theme.sunShapeGraphLineWidth = value
                            from: 1
                            to: 10
                        }
                        ColorPickerField {
                            label: "Color"
                            Layout.fillWidth: true
                            color: App.theme.sunShapeGraphLineColor
                            onColorChanged: App.theme.sunShapeGraphLineColor = color
                        }
                    }
                }
            }
        }

        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: parent.height - graphHeader.height - parent.spacing

            GraphsView {
                id: graphView
                anchors.fill: parent
                marginTop: 10
                marginLeft: 25
                marginBottom: 12

                theme: GraphsTheme {
                    colorScheme: GraphsTheme.ColorScheme.Dark
                    backgroundVisible: false
                    plotAreaBackgroundVisible: false
                    seriesColors: [App.theme.sunShapeGraphLineColor]
                    borderColors: [App.theme.fontColor, App.theme.fontColor]
                    grid.mainColor: Qt.rgba(App.theme.fontColor.r, App.theme.fontColor.g, App.theme.fontColor.b, 0.6)
                    grid.subColor: Qt.rgba(App.theme.fontColor.r, App.theme.fontColor.g, App.theme.fontColor.b, 0.6)
                    grid.mainWidth: 1
                    grid.subWidth: 1
                    axisX.mainColor:  App.theme.fontColor
                    axisX.subColor:  Qt.rgba(App.theme.fontColor.r, App.theme.fontColor.g, App.theme.fontColor.b, 0.6)
                    axisX.labelTextColor: App.theme.fontColor
                    axisX.mainWidth: 1
                    axisX.subWidth: 1
                    axisY.mainColor: App.theme.fontColor
                    axisY.subColor:  Qt.rgba(App.theme.fontColor.r, App.theme.fontColor.g, App.theme.fontColor.b, 0.6)
                    axisY.labelTextColor: App.theme.fontColor
                    axisY.mainWidth: 1
                    axisY.subWidth: 1
                    labelTextColor: App.theme.fontColor
                }

                axisX: ValueAxis {
                    id: xAxis
                    titleVisible: false
                    min: root.model ? root.model.x_axis_from : root.xMin
                    max: root.model ? root.model.x_axis_to : root.xMax
                }

                axisY: ValueAxis {
                    id: yAxis
                    titleVisible: false
                    min: root.model ? root.model.y_axis_from : root.yMin
                    max: root.model ? root.model.y_axis_to : root.yMax
                }

                LineSeries {
                    id: lineSeries
                    width: App.theme.sunShapeGraphLineWidth
                    name: ""
                }
            }

            Label {
                id: yAxisLabel
                text: root.yAxisTitle
                rotation: -90
                anchors.left: parent.left
                anchors.leftMargin: 0
                anchors.verticalCenter: graphView.verticalCenter
                anchors.verticalCenterOffset: -(xAxisLabel.height + xAxisLabel.anchors.bottomMargin)
            }

            Label {
                id: xAxisLabel
                text: root.xAxisTitle
                anchors.bottom: parent.bottom
                anchors.bottomMargin: 0
                anchors.horizontalCenter: graphView.horizontalCenter
                anchors.horizontalCenterOffset: yAxisLabel.height + yAxisLabel.anchors.leftMargin
            }
        }
    }
}
