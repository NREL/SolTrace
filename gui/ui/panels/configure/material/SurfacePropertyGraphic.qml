import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Shapes

import SolTrace

Item {
    id: root

    property real reflectance: 0.85   // ρ
    property real transmittance: 0.10 // τ
    property real nFront: 1.0
    property real nBack: 1.5
    property real slopeErrorMrad: 8.0
    property real specularityErrorMrad: 3.0

    property real incidentRayAngle: 33

    property bool exaggerateAngle: true

    property real radius: Math.min(root.height * 0.5, root.width * .4)

    property real band_radius: radius * .5

    Rectangle {
        anchors.fill: parent
        radius: 10
        color: Material.backgroundColor
    }

    Rectangle {
        anchors.top: parent.top
        anchors.left: parent.left
        anchors.right: parent.right

        anchors.bottom: parent.verticalCenter

        color: Qt.alpha(Material.color(Material.Blue), .1)

        topLeftRadius: 10
        topRightRadius: 10
    }

    Rectangle {
        anchors.bottom: parent.bottom
        anchors.left: parent.left
        anchors.right: parent.right

        anchors.top: parent.verticalCenter

        color: Qt.alpha(Material.color(Material.Green), .1)

        bottomLeftRadius: 10
        bottomRightRadius: 10
    }

    Rectangle {
        id: surface_level
        color: Material.foreground

        opacity: .75

        anchors.centerIn: parent
        width: parent.width * .75
        height: 3

        Label {
            text: "Surface"
            anchors.bottom: parent.top
            anchors.right: parent.right
            opacity: .5
            font.pointSize: 10
        }
    }

    component Ray : Rectangle {
        id: ray
        antialiasing: true

        property string title

        property real ray_angle: 0.0

        property bool other_side: false

        property bool label_lower: false

        width: root.radius
        height: 2

        anchors.right: other_side ? undefined : parent.horizontalCenter
        anchors.left:  other_side ? parent.horizontalCenter : undefined
        anchors.verticalCenter: parent.verticalCenter

        transform: [
            Rotation {
                origin.x: ray.other_side ? 0 : ray.width
                origin.y: ray.height / 2
                angle: ray_angle
            }
        ]

        Label {
            text: ray.title
            anchors.left: parent.other_side ? undefined : parent.left
            anchors.right: parent.other_side ? parent.right : undefined
            anchors.bottom: label_lower ? undefined : parent.top
            //anchors.bottomMargin: 4
            anchors.top: label_lower ? parent.bottom : undefined
            anchors.margins: 5

            //anchors.topMargin:
            font.pointSize: 11
            color: ray.color
        }

        Rectangle {
            id: arrow_1
            height: parent.height
            width: 5

            antialiasing: true

            anchors.right: parent.right
            anchors.verticalCenter: parent.verticalCenter

            color: ray.color

            transform: [
                Rotation {
                    origin.x: arrow_1.width
                    origin.y: arrow_1.height / 2
                    angle: 45
                }
            ]
        }

        Rectangle {
            id: arrow_2
            height: parent.height
            width: 5

            antialiasing: true

            anchors.right: parent.right
            anchors.verticalCenter: parent.verticalCenter

            color: ray.color

            transform: [
                Rotation {
                    origin.x: arrow_2.width
                    origin.y: arrow_2.height / 2
                    angle: -45
                }
            ]
        }
    }


    Ray {
        id: incident_ray

        title: "Incident"

        color: Material.color(Material.Orange)

        ray_angle: root.incidentRayAngle
    }

    Ray {
        id: reflect_ray

        title: "ρ"

        other_side: true

        color: Material.color(Material.Blue)

        ray_angle: -root.incidentRayAngle
    }

    Ray {
        id: transmit_ray

        title: "τ"

        other_side: true

        label_lower: true

        color: Material.color(Material.Green)

        ray_angle: root.incidentRayAngle
    }

    // We use two lines per band to avoid weird dash artifacts

    Shape {
        id: error_shape

        anchors.centerIn: parent
        width: root.radius * 2
        height: root.radius * 2

        preferredRendererType: Shape.CurveRenderer

        property vector2d center_pin: Qt.vector2d(width / 2, height / 2)

        // Angle convention:
        // 0 deg = right, -90 deg = up, +90 deg = down
        property real normalAngleDeg: -90
        property real reflectedAngleDeg: -root.incidentRayAngle

        property real slopeHalfAngleDeg: root.slopeErrorMrad * 180 / Math.PI / 1000 * (root.exaggerateAngle ? 10 : 1)
        property real specularHalfAngleDeg: root.specularityErrorMrad * 180 / Math.PI / 1000 * (root.exaggerateAngle ? 10 : 1)

        function degToRad(deg) {
            return deg * Math.PI / 180
        }

        function point_x(angleDeg) {
            return Math.cos(degToRad(angleDeg)) * root.band_radius + center_pin.x
        }

        function point_y(angleDeg) {
            return Math.sin(degToRad(angleDeg)) * root.band_radius + center_pin.y
        }

        ShapePath {
            strokeWidth: 2
            strokeColor: "grey"
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 2]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: Math.cos(error_shape.degToRad(-90)) * root.radius * .75 + error_shape.center_pin.x
                y: Math.sin(error_shape.degToRad(-90)) * root.radius * .75 + error_shape.center_pin.y
            }
        }

        // ----------------------------
        // Yellow slope-error band around vertical normal
        // ----------------------------

        ShapePath {
            id: surface_error_wedge

            strokeWidth: 0
            strokeColor: "transparent"
            fillColor: Qt.alpha(Material.color(Material.Yellow), 0.1)

            readonly property real a0: error_shape.normalAngleDeg
                                       - error_shape.slopeHalfAngleDeg
            readonly property real a1: error_shape.normalAngleDeg
                                       + error_shape.slopeHalfAngleDeg

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: error_shape.point_x(surface_error_wedge.a0)
                y: error_shape.point_y(surface_error_wedge.a0)
            }

            PathArc {
                x: error_shape.point_x(surface_error_wedge.a1)
                y: error_shape.point_y(surface_error_wedge.a1)

                radiusX: root.band_radius
                radiusY: root.band_radius

                useLargeArc: error_shape.slopeHalfAngleDeg > 90
            }

            PathLine {
                x: error_shape.center_pin.x
                y: error_shape.center_pin.y
            }
        }

        ShapePath {
            strokeWidth: 2
            strokeColor: Material.color(Material.Yellow)
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 2]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                id: slope_error_left
                x: error_shape.point_x(
                       error_shape.normalAngleDeg - error_shape.slopeHalfAngleDeg,
                       )
                y: error_shape.point_y(
                       error_shape.normalAngleDeg - error_shape.slopeHalfAngleDeg,
                       )
            }
        }

        ShapePath {
            strokeWidth: 2
            strokeColor: Material.color(Material.Yellow)
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 2]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: error_shape.point_x(
                       error_shape.normalAngleDeg + error_shape.slopeHalfAngleDeg,
                       )
                y: error_shape.point_y(
                       error_shape.normalAngleDeg + error_shape.slopeHalfAngleDeg,
                       )
            }
        }

        // Optional center line for the vertical normal
        ShapePath {
            strokeWidth: 1
            strokeColor: Qt.alpha(Material.color(Material.Yellow), 0.5)
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 3]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: error_shape.point_x(error_shape.normalAngleDeg)
                y: error_shape.point_y(error_shape.normalAngleDeg)
            }
        }



        // ----------------------------
        // Pink specular-error band around reflected ray
        // ----------------------------

        ShapePath {
            id: specular_error_wedge

            strokeWidth: 0
            strokeColor: "transparent"
            fillColor: Qt.alpha(Material.color(Material.Pink), 0.1)

            readonly property real a0: error_shape.reflectedAngleDeg
                                       - error_shape.specularHalfAngleDeg
            readonly property real a1: error_shape.reflectedAngleDeg
                                       + error_shape.specularHalfAngleDeg

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: error_shape.point_x(specular_error_wedge.a0)
                y: error_shape.point_y(specular_error_wedge.a0)
            }

            PathArc {
                x: error_shape.point_x(specular_error_wedge.a1)
                y: error_shape.point_y(specular_error_wedge.a1)

                radiusX: root.band_radius
                radiusY: root.band_radius

                useLargeArc: error_shape.specularHalfAngleDeg > 90
            }

            PathLine {
                x: error_shape.center_pin.x
                y: error_shape.center_pin.y
            }
        }

        ShapePath {
            strokeWidth: 2
            strokeColor: Material.color(Material.Pink)
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 2]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                x: error_shape.point_x(
                       error_shape.reflectedAngleDeg - error_shape.specularHalfAngleDeg,
                       )
                y: error_shape.point_y(
                       error_shape.reflectedAngleDeg - error_shape.specularHalfAngleDeg,
                       )
            }
        }

        ShapePath {
            strokeWidth: 2
            strokeColor: Material.color(Material.Pink)
            strokeStyle: ShapePath.DashLine
            fillColor: "transparent"
            dashPattern: [1, 2]

            startX: error_shape.center_pin.x
            startY: error_shape.center_pin.y

            PathLine {
                id: spec_error_left

                x: error_shape.point_x(
                       error_shape.reflectedAngleDeg + error_shape.specularHalfAngleDeg,
                       )
                y: error_shape.point_y(
                       error_shape.reflectedAngleDeg + error_shape.specularHalfAngleDeg,
                       )
            }
        }
    }

    property vector2d slope_error_placement: mapFromItem(error_shape,
                                             Qt.vector2d(slope_error_left.x, slope_error_left.y))

    Label {
        color: Material.color(Material.Yellow)
        text : "σ<sub>slope</sub>"
        textFormat: Label.RichText

        x: parent.slope_error_placement.x - width
        y: parent.slope_error_placement.y - 20
    }

    property vector2d spec_error_placement: mapFromItem(error_shape,
                                             Qt.vector2d(spec_error_left.x, spec_error_left.y))

    Label {
        color: Material.color(Material.Pink)
        text : "σ<sub>spec</sub>"
        textFormat: Label.RichText

        x: parent.spec_error_placement.x
        y: parent.spec_error_placement.y
    }

    Label {
        color: Material.color(Material.Grey)
        text: "<em>n</em> Front = " + root.nFront
        textFormat: Label.RichText

        x: parent.width * .05
        anchors.bottom: parent.verticalCenter
        anchors.bottomMargin: 3
    }

    Label {
        color: Material.color(Material.Grey)
        text: "<em>n</em> Back = " + root.nBack
        textFormat: Label.RichText

        x: parent.width * .05
        anchors.top: parent.verticalCenter
        anchors.topMargin: 3
    }

    Label {
        color: Material.color(Material.Green)
        text: root.transmittance
        textFormat: Label.RichText

        anchors.right: parent.right
        anchors.rightMargin: 10
        y: parent.height * .75 - height / 2
    }

    Label {
        color: Material.color(Material.Blue)
        text: root.reflectance
        textFormat: Label.RichText

        anchors.right: parent.right
        anchors.rightMargin: 10
        y: parent.height * .25 - height / 2
    }

    STIconButton {
        text: "\uf013"

        onClicked: graphic_pop.open()

        anchors.bottom: parent.bottom
        anchors.left: parent.left

        anchors.margins: 3

        STPopup {
            id: graphic_pop

            ColumnLayout {
                anchors.fill: parent
                STSwitch {
                    id: ex_angles
                    text: "Exaggerate angles"
                    checked: root.exaggerateAngle

                    onToggled: {
                        root.exaggerateAngle = !root.exaggerateAngle
                    }
                }
            }
        }
    }
}
