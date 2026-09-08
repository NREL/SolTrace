import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Layouts

Item {
    id: root
    width: 980
    height: 620

    readonly property var surfaceNames: [
        "CONE",
        "CYLINDER",
        "FLAT",
        "PARABOLA",
        "SPHERE",
        "HYPER",
        "GENERAL_SPENCER_MURTY",
        "TORUS"
    ]

    GridLayout {
        anchors.fill: parent
        anchors.margins: 16
        columns: 4
        rowSpacing: 12
        columnSpacing: 12

        Repeater {
            model: root.surfaceNames

            delegate: Frame {
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.preferredWidth: 220
                Layout.preferredHeight: 180
                padding: 8

                // background: Rectangle {
                //     radius: 8
                //     color: "#1c2430"
                //     border.color: "#3a4658"
                //     border.width: 1
                // }

                ColumnLayout {
                    anchors.fill: parent
                    spacing: 6

                    Label {
                        text: modelData
                        color: "#dfe7f2"
                        font.bold: true
                        font.pointSize: 14
                        Layout.alignment: Qt.AlignHCenter
                    }

                    Canvas {
                        id: canvas
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        antialiasing: true

                        onPaint: {
                            const ctx = getContext("2d");
                            const w = width;
                            const h = height;

                            ctx.reset();
                            ctx.clearRect(0, 0, w, h);

                            // background
                            ctx.fillStyle = "#11161d";
                            roundRect(ctx, 0, 0, w, h, 8);
                            ctx.fill();

                            drawGrid(ctx, w, h);
                            drawAxes(ctx, w, h);

                            switch (modelData) {
                            case "CONE":
                                drawCone(ctx, w, h);
                                break;
                            case "CYLINDER":
                                drawCylinder(ctx, w, h);
                                break;
                            case "FLAT":
                                drawFlat(ctx, w, h);
                                break;
                            case "PARABOLA":
                                drawParabola(ctx, w, h);
                                break;
                            case "SPHERE":
                                drawSphere(ctx, w, h);
                                break;
                            case "HYPER":
                                drawHyperbola(ctx, w, h);
                                break;
                            case "GENERAL_SPENCER_MURTY":
                                drawGeneralSpencerMurty(ctx, w, h);
                                break;
                            case "TORUS":
                                drawTorus(ctx, w, h);
                                break;
                            }
                        }

                        Component.onCompleted: requestPaint()
                        onWidthChanged: requestPaint()
                        onHeightChanged: requestPaint()

                        function roundRect(ctx, x, y, w, h, r) {
                            ctx.beginPath();
                            ctx.moveTo(x + r, y);
                            ctx.lineTo(x + w - r, y);
                            ctx.quadraticCurveTo(x + w, y, x + w, y + r);
                            ctx.lineTo(x + w, y + h - r);
                            ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
                            ctx.lineTo(x + r, y + h);
                            ctx.quadraticCurveTo(x, y + h, x, y + h - r);
                            ctx.lineTo(x, y + r);
                            ctx.quadraticCurveTo(x, y, x + r, y);
                            ctx.closePath();
                        }

                        function drawGrid(ctx, w, h) {
                            ctx.save();
                            ctx.strokeStyle = "#1f2a36";
                            ctx.lineWidth = 1;

                            const step = 20;
                            for (let x = step; x < w; x += step) {
                                ctx.beginPath();
                                ctx.moveTo(x, 0);
                                ctx.lineTo(x, h);
                                ctx.stroke();
                            }
                            for (let y = step; y < h; y += step) {
                                ctx.beginPath();
                                ctx.moveTo(0, y);
                                ctx.lineTo(w, y);
                                ctx.stroke();
                            }
                            ctx.restore();
                        }

                        function drawAxes(ctx, w, h) {
                            ctx.save();
                            ctx.strokeStyle = "#516176";
                            ctx.lineWidth = 1.2;

                            const cx = w * 0.5;
                            const cy = h * 0.56;

                            ctx.beginPath();
                            ctx.moveTo(10, cy);
                            ctx.lineTo(w - 10, cy);
                            ctx.stroke();

                            ctx.beginPath();
                            ctx.moveTo(cx, 10);
                            ctx.lineTo(cx, h - 10);
                            ctx.stroke();

                            ctx.restore();
                        }

                        function strokeProfile(ctx) {
                            ctx.lineWidth = 3;
                            ctx.strokeStyle = "#7fd3ff";
                            ctx.stroke();

                            ctx.lineWidth = 1;
                            ctx.strokeStyle = "#c7f0ff";
                            ctx.stroke();
                        }

                        function drawCone(ctx, w, h) {
                            const cx = w * 0.5;
                            const topY = h * 0.22;
                            const baseY = h * 0.82;
                            const halfBase = w * 0.22;

                            ctx.beginPath();
                            ctx.moveTo(cx, topY);
                            ctx.lineTo(cx - halfBase, baseY);
                            ctx.lineTo(cx + halfBase, baseY);
                            ctx.closePath();
                            strokeProfile(ctx);

                            ctx.beginPath();
                            ctx.ellipse(cx, baseY, halfBase, 10, 0, 0, Math.PI * 2);
                            ctx.lineWidth = 2;
                            ctx.strokeStyle = "#7fd3ff";
                            ctx.stroke();
                        }

                        function drawCylinder(ctx, w, h) {
                            const cx = w * 0.5;
                            const topY = h * 0.22;
                            const botY = h * 0.82;
                            const rx = w * 0.18;

                            ctx.beginPath();
                            ctx.moveTo(cx - rx, topY);
                            ctx.lineTo(cx - rx, botY);
                            ctx.moveTo(cx + rx, topY);
                            ctx.lineTo(cx + rx, botY);
                            strokeProfile(ctx);

                            ctx.beginPath();
                            ctx.ellipse(cx, topY, rx, 10, 0, 0, Math.PI * 2);
                            ctx.strokeStyle = "#7fd3ff";
                            ctx.lineWidth = 2;
                            ctx.stroke();

                            ctx.beginPath();
                            ctx.ellipse(cx, botY, rx, 10, 0, 0, Math.PI * 2);
                            ctx.strokeStyle = "#7fd3ff";
                            ctx.lineWidth = 2;
                            ctx.stroke();
                        }

                        function drawFlat(ctx, w, h) {
                            const y = h * 0.56;

                            ctx.beginPath();
                            ctx.moveTo(w * 0.18, y);
                            ctx.lineTo(w * 0.82, y);
                            strokeProfile(ctx);
                        }

                        function drawParabola(ctx, w, h) {
                            const cx = w * 0.5;
                            const topY = h * 0.2;
                            const bottomY = h * 0.82;
                            const maxHalfWidth = w * 0.28;

                            ctx.beginPath();
                            for (let i = 0; i <= 100; ++i) {
                                const t = i / 100.0;      // 0..1 down the shape
                                const y = topY + t * (bottomY - topY);
                                const x = maxHalfWidth * Math.sqrt(t);
                                const pxL = cx - x;
                                const pxR = cx + x;

                                if (i === 0)
                                    ctx.moveTo(pxL, y);
                                else
                                    ctx.lineTo(pxL, y);
                            }
                            for (let i = 100; i >= 0; --i) {
                                const t = i / 100.0;
                                const y = topY + t * (bottomY - topY);
                                const x = maxHalfWidth * Math.sqrt(t);
                                ctx.lineTo(cx + x, y);
                            }
                            ctx.closePath();
                            strokeProfile(ctx);

                            // focus marker
                            ctx.beginPath();
                            ctx.arc(cx, h * 0.37, 3, 0, Math.PI * 2);
                            ctx.fillStyle = "#ffd166";
                            ctx.fill();
                        }

                        function drawSphere(ctx, w, h) {
                            const cx = w * 0.5;
                            const cy = h * 0.56;
                            const r = Math.min(w, h) * 0.28;

                            ctx.beginPath();
                            ctx.arc(cx, cy, r, 0, Math.PI * 2);
                            strokeProfile(ctx);

                            // equator hint
                            ctx.beginPath();
                            ctx.ellipse(cx, cy, r, r * 0.28, 0, 0, Math.PI * 2);
                            ctx.lineWidth = 1.5;
                            ctx.strokeStyle = "#95e1ff";
                            ctx.stroke();
                        }

                        function drawHyperbola(ctx, w, h) {
                            const cx = w * 0.5;
                            const cy = h * 0.56;
                            const a = w * 0.12;
                            const b = h * 0.16;

                            ctx.beginPath();
                            for (let t = -1.15; t <= 1.15; t += 0.02) {
                                const x = a * Math.cosh(t);
                                const y = b * Math.sinh(t);
                                const px = cx - x;
                                const py = cy + y;
                                if (t === -1.15) ctx.moveTo(px, py);
                                else ctx.lineTo(px, py);
                            }
                            for (let t = 1.15; t >= -1.15; t -= 0.02) {
                                const x = a * Math.cosh(t);
                                const y = b * Math.sinh(t);
                                ctx.lineTo(cx + x, cy + y);
                            }
                            ctx.closePath();
                            strokeProfile(ctx);
                        }

                        function drawGeneralSpencerMurty(ctx, w, h) {
                            // Generic freeform / aspheric optical profile sketch
                            // for a "general" surface tile.
                            const x0 = w * 0.18;
                            const x1 = w * 0.82;

                            ctx.beginPath();
                            for (let i = 0; i <= 120; ++i) {
                                const t = i / 120.0;
                                const x = x0 + t * (x1 - x0);

                                // smooth freeform-ish profile
                                const y = h * (0.72
                                             - 0.30 * t
                                             + 0.10 * t * t
                                             - 0.22 * Math.pow(t - 0.58, 4)
                                             + 0.03 * Math.sin(t * Math.PI * 2.4));

                                if (i === 0) ctx.moveTo(x, y);
                                else ctx.lineTo(x, y);
                            }
                            strokeProfile(ctx);

                            // receiver / target marker
                            ctx.beginPath();
                            ctx.arc(w * 0.74, h * 0.26, 4, 0, Math.PI * 2);
                            ctx.fillStyle = "#ffd166";
                            ctx.fill();
                        }

                        function drawTorus(ctx, w, h) {
                            const cx = w * 0.5;
                            const cy = h * 0.56;
                            const outerRx = w * 0.28;
                            const outerRy = h * 0.22;
                            const innerRx = w * 0.13;
                            const innerRy = h * 0.09;

                            ctx.beginPath();
                            ctx.ellipse(cx, cy, outerRx, outerRy, 0, 0, Math.PI * 2);
                            strokeProfile(ctx);

                            ctx.beginPath();
                            ctx.ellipse(cx, cy, innerRx, innerRy, 0, 0, Math.PI * 2);
                            ctx.lineWidth = 2.2;
                            ctx.strokeStyle = "#7fd3ff";
                            ctx.stroke();
                        }
                    }
                }
            }
        }
    }
}
