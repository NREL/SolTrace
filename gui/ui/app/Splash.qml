import QtQuick

// Intentionally minimal to speed up creation
Window {
    id: splash_window

    title: "SolTrace"

    width: 400
    height: 256
    visible: true

    color: "transparent"

    modality: Qt.ApplicationModal
    flags: Qt.SplashScreen | Qt.FramelessWindowHint

    x: (Screen.width - width) / 2
    y: (Screen.height - height) / 2

    Rectangle {
        anchors.fill: parent
        gradient: Gradient {
            GradientStop { position: 0.0; color: "#aaaaaa" }
            GradientStop { position: 1.0; color: "#40525e" }
        }
        //color: "#40525e"
        radius: 10

        Image {
            anchors.fill: parent
            anchors.margins: 12
            fillMode: Image.PreserveAspectFit
            source: "qrc:/assets/images/CSP-SolTrace Logo-FY26-DB-Final_SolTrace_Logo_Horizontal_White.png"
            mipmap: true
        }

        // BusyIndicator {
        //     anchors.right: parent.right
        //     anchors.bottom: parent.bottom
        //     width: 24
        //     height: 24
        // }
    }
}
