import QtQuick
import QtQuick.Shapes
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material

import SolTrace

Item {
    id: root

    opacity: .75

    property var model: App.sun.shape.current_distribution
    property real max_angle: 1.0
    property real max_intensity: 0.0
    property real half_max_angle: 0.0
    property bool has_distribution: false
    property var current_gradient: null
    readonly property int max_gradient_stops: 100

    function computeMaxAngle() {
        if (!model) {
            return 1.0
        }

        var maxAngle = 0.0

        for (let i = 0; i < model.rowCount(); i++) {
            const angle = Math.max(0.0, model.data(model.index(i, 0)))

            maxAngle = Math.max(maxAngle, angle)
        }

        return Math.max(maxAngle, 1.0e-9)
    }

    property var color_map: [
        [0.101441, 0.200110, 0.700194],
        [0.125813, 0.244322, 0.678320],
        [0.147947, 0.288802, 0.655913],
        [0.172248, 0.337633, 0.631277],
        [0.202374, 0.390426, 0.601219],
        [0.240659, 0.436433, 0.563176],
        [0.283755, 0.480545, 0.523580],
        [0.335316, 0.532757, 0.495080],
        [0.394629, 0.592751, 0.477439],
        [0.454835, 0.653426, 0.462265],
        [0.517878, 0.716845, 0.446543],
        [0.584205, 0.783089, 0.430163],
        [0.661151, 0.851901, 0.414009],
        [0.754686, 0.908919, 0.404012],
        [0.876163, 0.957911, 0.400588],
        [1.000000, 0.999989, 0.400094],
    ]

    function sample_color_map(t) {
        // Clamp normalized input to [0, 1]
        t = Math.max(0, Math.min(1, t));

        const n = color_map.length - 1;
        const x = t * n;
        const i = Math.floor(x);
        const f = x - i;

        if (i >= n) return color_map[n];

        const a = color_map[i];
        const b = color_map[i + 1];

        return [
                    a[0] + (b[0] - a[0]) * f,
                    a[1] + (b[1] - a[1]) * f,
                    a[2] + (b[2] - a[2]) * f,
                ];
    }

    function sample_color_mapCss(t) {
        var [r, g, b] = sample_color_map(t);
        r = Math.round(r * 255)
        g = Math.round(g * 255)
        b = Math.round(b * 255)
        return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
    }

    function formatAngle(angle) {
        const value = Math.abs(angle)
        if (value >= 10.0) {
            return value.toFixed(1) + " mrad"
        }
        if (value >= 1.0) {
            return value.toFixed(2) + " mrad"
        }
        return value.toPrecision(2) + " mrad"
    }

    function computeHalfMaxAngle(samples, peakIntensity) {
        if (samples.length === 0 || peakIntensity <= 0.0) {
            return 0.0
        }

        const halfIntensity = peakIntensity * 0.5
        var previous = samples[0]

        if (previous.intensity <= halfIntensity) {
            return previous.angle
        }

        for (let i = 1; i < samples.length; ++i) {
            const current = samples[i]
            if (current.intensity > halfIntensity) {
                previous = current
                continue
            }

            const delta = current.intensity - previous.intensity
            if (Math.abs(delta) <= 1.0e-12) {
                return current.angle
            }

            const t = (halfIntensity - previous.intensity) / delta
            return previous.angle + (current.angle - previous.angle) * t
        }

        return samples[samples.length - 1].angle
    }

    function normalizedAngle(angle) {
        if (max_angle <= 0.0) {
            return 0.0
        }
        return Math.max(0.0, Math.min(1.0, Math.abs(angle) / max_angle))
    }

    function angleX(angle) {
        return chart.width / 2 + normalizedAngle(angle) * shape_path.diskRadius
    }

    function destroyRadialGradient(gradient) {
        if (!gradient) {
            return
        }

        const oldStops = gradient.stops
        gradient.stops = []
        for (let i = 0; i < oldStops.length; ++i) {
            oldStops[i].destroy()
        }
        gradient.destroy()
    }

    function addGradientStop(stops, position, color) {
        const stop = Qt.createQmlObject(
                       'import QtQuick; GradientStop {}',
                       root,
                       "sunPreviewGradientStop")
        stop.position = Math.max(0.0, Math.min(1.0, position))
        stop.color = color
        stops.push(stop)
    }

    function createRadialGradient(stops) {
        const gradient = Qt.createQmlObject(
                         'import QtQuick.Shapes; RadialGradient {}',
                         root,
                         "sunPreviewRadialGradient")
        gradient.centerX = Qt.binding(function() { return chart.width / 2 })
        gradient.centerY = Qt.binding(function() { return chart.height / 2 })
        gradient.focalX = Qt.binding(function() { return gradient.centerX })
        gradient.focalY = Qt.binding(function() { return gradient.centerY })
        gradient.centerRadius = Qt.binding(function() { return shape_path.diskRadius })
        gradient.stops = stops
        return gradient
    }

    function replaceRadialGradient(stops) {
        const oldGradient = current_gradient
        const gradient = createRadialGradient(stops)
        current_gradient = gradient
        shape_path.fillGradient = gradient
        destroyRadialGradient(oldGradient)
    }

    function sampledGradientSamples(samples) {
        if (samples.length <= max_gradient_stops) {
            return samples
        }

        const sampled = []
        const scale = (samples.length - 1) / (max_gradient_stops - 1)
        for (let i = 0; i < max_gradient_stops; ++i) {
            sampled.push(samples[Math.round(i * scale)])
        }
        return sampled
    }

    function rebuildGradientStops() {
        max_angle = computeMaxAngle()

        const stops = []
        if (!model || model.rowCount() === 0) {
            has_distribution = false
            max_intensity = 0.0
            half_max_angle = 0.0
            addGradientStop(stops, 0.0, sample_color_mapCss(0.0))
            addGradientStop(stops, 1.0, sample_color_mapCss(0.0))
            replaceRadialGradient(stops)
            return
        }

        const samples = []
        var maxIntensity = 0.0
        for (let i = 0; i < model.rowCount(); ++i) {
            const angle = Math.max(0.0, model.data(model.index(i, 0)))
            const intensity = Math.max(0.0, model.data(model.index(i, 1)))

            maxIntensity = Math.max(maxIntensity, intensity)
            samples.push({
                             "angle": angle,
                             "position": angle / max_angle,
                             "intensity": intensity
                         })
        }

        samples.sort(function(a, b) {
            return a.position - b.position
        })

        const mergedSamples = []
        for (let j = 0; j < samples.length; ++j) {
            const previous = mergedSamples[mergedSamples.length - 1]
            if (previous && Math.abs(previous.position - samples[j].position) < 1.0e-6) {
                previous.intensity = Math.max(previous.intensity, samples[j].intensity)
            } else {
                mergedSamples.push(samples[j])
            }
        }

        has_distribution = true
        max_intensity = maxIntensity
        half_max_angle = computeHalfMaxAngle(mergedSamples, maxIntensity)

        const scale = maxIntensity > 0.0 ? maxIntensity : 1.0
        const gradientSamples = sampledGradientSamples(mergedSamples)
        for (let k = 0; k < gradientSamples.length; ++k) {
            addGradientStop(stops,
                            gradientSamples[k].position,
                            sample_color_mapCss(gradientSamples[k].intensity / scale))
        }

        if (stops[0].position > 0.0) {
            addGradientStop(stops, 0.0, stops[0].color)
        }
        if (stops[stops.length - 1].position < 1.0) {
            addGradientStop(stops, 1.0, stops[stops.length - 1].color)
        }

        stops.sort(function(a, b) {
            return a.position - b.position
        })

        replaceRadialGradient(stops)
    }

    onModelChanged: rebuildGradientStops()

    Component.onCompleted: rebuildGradientStops()
    Component.onDestruction: destroyRadialGradient(current_gradient)

    Connections {
        target: root.model

        function onChanged() {
            root.rebuildGradientStops()
        }
    }

    Connections {
        target: App.sun.shape

        function onChanged() {
            root.rebuildGradientStops()
        }
    }

    RowLayout {
        anchors.fill: parent
        spacing: 10

        Item {
            id: sourceGraphic

            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.preferredWidth: 1

            property real graphicRadius: Math.max(0, Math.min(width, height) * 0.34)
            property real centerX: width / 2
            property real centerY: height / 2

            Repeater {
                model: 18

                Shape {
                    required property int index

                    anchors.fill: sourceGraphic
                    visible: App.sun.type === SunModule.PointSource
                    preferredRendererType: Shape.CurveRenderer

                    property real theta: 2.0 * Math.PI * index / 18.0
                    property real innerRadius: sourceGraphic.graphicRadius * 0.34
                    property real outerRadius: sourceGraphic.graphicRadius * 0.95

                    ShapePath {
                        strokeWidth: 5
                        strokeColor: Qt.alpha(App.theme.fontColor, 0.78)
                        fillColor: "transparent"
                        capStyle: ShapePath.RoundCap

                        startX: sourceGraphic.centerX + Math.cos(theta) * innerRadius
                        startY: sourceGraphic.centerY + Math.sin(theta) * innerRadius

                        PathLine {
                            x: sourceGraphic.centerX + Math.cos(theta) * outerRadius
                            y: sourceGraphic.centerY + Math.sin(theta) * outerRadius
                        }
                    }
                }
            }

            Shape {
                anchors.fill: parent
                visible: App.sun.type === SunModule.PointSource
                preferredRendererType: Shape.CurveRenderer

                ShapePath {
                    strokeWidth: 5
                    strokeColor: App.theme.fontColor
                    fillColor: Qt.alpha(App.theme.fontColor, 0.10)

                    startX: sourceGraphic.centerX + sourceGraphic.graphicRadius * 0.34
                    startY: sourceGraphic.centerY

                    PathAngleArc {
                        centerX: sourceGraphic.centerX
                        centerY: sourceGraphic.centerY
                        radiusX: sourceGraphic.graphicRadius * 0.34
                        radiusY: radiusX
                        startAngle: 0
                        sweepAngle: 360
                    }
                }
            }

            component DirectionArrow: Shape {
                property real yOffset: 0.0

                anchors.fill: sourceGraphic
                visible: App.sun.type === SunModule.Directional
                preferredRendererType: Shape.CurveRenderer

                property real startXPos: sourceGraphic.centerX - sourceGraphic.graphicRadius * 0.74
                property real endXPos: sourceGraphic.centerX + sourceGraphic.graphicRadius * 0.74
                property real yPos: sourceGraphic.centerY + yOffset
                property real headSize: Math.max(12, sourceGraphic.graphicRadius * 0.20)

                ShapePath {
                    strokeWidth: 6
                    strokeColor: App.theme.fontColor
                    fillColor: "transparent"
                    capStyle: ShapePath.RoundCap
                    joinStyle: ShapePath.RoundJoin

                    startX: startXPos
                    startY: yPos

                    PathLine {
                        x: endXPos
                        y: yPos
                    }
                }

                ShapePath {
                    strokeWidth: 6
                    strokeColor: App.theme.fontColor
                    fillColor: "transparent"
                    capStyle: ShapePath.RoundCap
                    joinStyle: ShapePath.RoundJoin

                    startX: endXPos - headSize
                    startY: yPos - headSize * 0.55

                    PathLine {
                        x: endXPos
                        y: yPos
                    }

                    PathLine {
                        x: endXPos - headSize
                        y: yPos + headSize * 0.55
                    }
                }
            }

            DirectionArrow {
                yOffset: -sourceGraphic.graphicRadius * 0.42
            }

            DirectionArrow {
                yOffset: 0
            }

            DirectionArrow {
                yOffset: sourceGraphic.graphicRadius * 0.42
            }
        }

        Item {
            id: chart

            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.preferredWidth: 1

            Shape {
                anchors.fill: parent

                preferredRendererType: Shape.CurveRenderer

                ShapePath {
                    id: shape_path
                    property real diskRadius: Math.max(0, Math.min(chart.width, chart.height) / 2 - 1)

                    startX: chart.width / 2 + diskRadius
                    startY: chart.height / 2
                    strokeWidth: 1
                    strokeColor: "white"
                    fillColor: "white"

                    PathAngleArc {
                        centerX: chart.width / 2
                        centerY: chart.height / 2
                        radiusX: shape_path.diskRadius
                        radiusY: radiusX
                        startAngle: 0
                        sweepAngle: 360
                    }
                }
            }

            component RadialTick: Shape {
                property real angle: 0.0

                anchors.fill: chart
                visible: root.has_distribution
                z: 1
                preferredRendererType: Shape.CurveRenderer

                property real tickRadius: root.normalizedAngle(angle) * shape_path.diskRadius

                ShapePath {
                    strokeWidth: 1
                    strokeColor: Qt.alpha(App.theme.fontColor, 0.46)
                    fillColor: "transparent"
                    strokeStyle: ShapePath.DashLine
                    dashPattern: [4, 4]

                    startX: chart.width / 2 + tickRadius
                    startY: chart.height / 2

                    PathAngleArc {
                        centerX: chart.width / 2
                        centerY: chart.height / 2
                        radiusX: tickRadius
                        radiusY: tickRadius
                        startAngle: 0
                        sweepAngle: 360
                    }
                }
            }

            RadialTick {
                angle: root.max_angle * 0.25
            }

            RadialTick {
                angle: root.max_angle * 0.5
            }

            RadialTick {
                angle: root.max_angle * 0.75
            }

            RadialTick {
                angle: root.max_angle
            }

            component AngleLabel: Item {
                property real angle: 0.0
                property real verticalOffset: 0.0

                visible: root.has_distribution
                z: 2

                width: Math.min(label.implicitWidth + 12, chart.width - 12)
                height: label.implicitHeight + 6
                x: Math.max(6, Math.min(chart.width - width - 6, root.angleX(angle) - width / 2))
                y: Math.max(6, Math.min(chart.height - height - 6, chart.height / 2 + verticalOffset - height / 2))

                Rectangle {
                    anchors.fill: parent
                    radius: 4
                    color: Qt.rgba(0, 0, 0, 0.42)
                    border.width: 1
                    border.color: Qt.alpha(App.theme.fontColor, 0.25)
                }

                Label {
                    id: label
                    anchors.fill: parent
                    anchors.leftMargin: 6
                    anchors.rightMargin: 6
                    color: App.theme.fontColor
                    font.pointSize: 10
                    horizontalAlignment: Text.AlignHCenter
                    verticalAlignment: Text.AlignVCenter
                    elide: Label.ElideRight
                    text: root.formatAngle(angle)
                }
            }

            AngleLabel {
                angle: root.max_angle
                verticalOffset: 14
            }
        }
    }
}
