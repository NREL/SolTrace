import QtQuick
import QtQuick.Controls.Material
import QtQuick.Layouts

import SolTrace

ScrollView {
    id: root
    Layout.fillWidth: true
    Layout.fillHeight: true
    contentWidth: availableWidth

    ColumnLayout {
        width: root.availableWidth
        spacing: 12

        ListModel {
            id: backendFeatureModel
            ListElement {
                name: "GPU Ray Tracing"
                icon: "\uf2db"
                description: "Massively parallel Monte Carlo ray tracing on GPU architecture. Full physical fidelity with order-of-magnitude speedups over legacy CPU-bound workflows."
            }
            ListElement {
                name: "Large-Scale Simulations"
                icon: "\uf0b2"
                description: "Robust handling of high-resolution, full-scale heliostat fields. Simulations that were previously impractical are now routine."
            }
            ListElement {
                name: "Open-Source & Extensible"
                icon: "\uf121"
                description: "Modular codebase with clear APIs. Customizable beyond proprietary software limits, lowering barriers to entry for researchers and developers."
            }
            ListElement {
                name: "SAM Integration"
                icon: "\uf0c1"
                description: "Coupling optical simulations with system-level techno-economic models like SAM for end-to-end CST performance and cost analysis."
            }
            ListElement {
                name: "Python API"
                icon: "\uf120"
                description: "Scriptable interface via pysoltrace for parametric studies, optimization loops, and integration with existing Python-based research pipelines."
            }
        }

        ListModel {
            id: guiFeatureModel
            ListElement {
                name: "3D Viewport"
                icon: "\uf1b2"
                description: "Interactive 3D visualization of optical geometries, ray intersections, and flux maps. Explore your system spatially instead of through tables."
            }
            ListElement {
                name: "Geometry Editor"
                icon: "\uf5ee"
                description: "Direct manipulation of stages, elements, and optical properties. Eliminates trial-and-error geometry generation from legacy workflows."
            }
            ListElement {
                name: "Flux Visualization"
                icon: "\uf06d"
                description: "Real-time flux map rendering and scatter plots. Quickly assess optical performance without exporting data to external tools."
            }
            ListElement {
                name: "Theming & Accessibility"
                icon: "\uf53f"
                description: "Customizable glass UI with adjustable colors, font sizes, and zoom levels. Adapts to user preferences and display environments."
            }
            ListElement {
                name: "Documentation"
                icon: "\uf02d"
                description: "Built-in multilingual documentation and help resources. Minimizes learning curve and provides guidance without leaving the application."
            }
        }

        Header {
            text: "New in SolTrace"
        }

        STPropertyPanel {
            title: "Backend Features"
            Layout.fillWidth: true
            CardGallery {
                model: backendFeatureModel
                Layout.fillWidth: true
                Layout.columnSpan: 2
                delegate: ColumnLayout {
                    property string name
                    property string icon
                    property string description
                    spacing: 5
                    RowLayout {
                        spacing: 8
                        Label {
                            text: icon
                            font.family: "Font Awesome 7 Free"
                            font.pointSize: 16
                        }
                        SubHeader { text: name }
                    }
                    Label {
                        Layout.fillWidth: true
                        text: description
                        wrapMode: Text.WordWrap
                    }
                }
            }
        }

        STPropertyPanel {
            title: "GUI Features"
            Layout.fillWidth: true
            CardGallery {
                model: guiFeatureModel
                Layout.fillWidth: true
                Layout.columnSpan: 2
                delegate: ColumnLayout {
                    property string name
                    property string icon
                    property string description
                    spacing: 5
                    RowLayout {
                        spacing: 8
                        Label {
                            text: icon
                            font.family: "Font Awesome 7 Free"
                            font.pointSize: 16
                        }
                        SubHeader { text: name }
                    }
                    Label {
                        Layout.fillWidth: true
                        text: description
                        wrapMode: Text.WordWrap
                    }
                }
            }
        }
    }
}
