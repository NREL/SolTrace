import QtQuick
import QtCore
import QtQuick3D
import QtQuick3D.Helpers
import QtQuick3D.AssetUtils

QtObject {
    id: theme

    property color defaultGlassColor: Qt.rgba(0, 0, 0, 0.15)
    property color defaultFontColor: "white"

    property color glassColor: defaultGlassColor
    readonly property color destructiveGlassColor: Qt.rgba(glassColor.r + 0.1, glassColor.g, glassColor.b, glassColor.a)
    readonly property color shadedGlassColor: theme.glassColorA(glassColor.a + 0.05)
    readonly property color dividerColor: Qt.rgba(1, 1, 1, 0.2)

    property color fontColor: defaultFontColor

    readonly property int defaultHeaderSize: 17
    readonly property int defaultSubHeaderSize: 16
    readonly property int defaultLabelSize: 13
    readonly property int defaultComboBarTextSize: 13
    readonly property int defaultPropertyPanelHeaderSize: 14
    readonly property int defaultNormalSize: 15

    property int _headerSize: defaultHeaderSize
    property int _subHeaderSize: defaultSubHeaderSize
    property int _labelSize: defaultLabelSize
    property int _comboBarTextSize: defaultComboBarTextSize
    property int _propertyPanelHeaderSize: defaultPropertyPanelHeaderSize
    property int _normalSize: defaultNormalSize

    property real zoomLevel: 1

    function calculateSize(baseSize) {
        return Math.max(1, Math.min(48, (zoomLevel - 1) * 4 + baseSize))
    }

    readonly property int headerSize: calculateSize(_headerSize)
    readonly property int subHeaderSize: calculateSize(_subHeaderSize)
    readonly property int labelSize: calculateSize(_labelSize)
    readonly property int comboBarTextSize: calculateSize(_comboBarTextSize)
    readonly property int propertyPanelHeaderSize: calculateSize(_propertyPanelHeaderSize)
    readonly property int normalSize: calculateSize(_normalSize)

    property color sunShapeGraphLineColor: "#E0D080"
    property int sunShapeGraphLineWidth: 3

    readonly property Settings settings: Settings {
        property alias glassColor: theme.glassColor
        property alias fontColor: theme.fontColor
        property alias zoomLevel: theme.zoomLevel
        property alias headerSize: theme._headerSize
        property alias subHeaderSize: theme._subHeaderSize
        property alias labelSize: theme._labelSize
        property alias comboBarTextSize: theme._comboBarTextSize
        property alias propertyPanelHeaderSize: theme._propertyPanelHeaderSize
        property alias normalSize: theme._normalSize
        property alias sunShapeGraphLineColor: theme.sunShapeGraphLineColor
        property alias sunShapeGraphLineWidth: theme.sunShapeGraphLineWidth
    }

    function glassColorA(alpha) {
        return Qt.rgba(glassColor.r, glassColor.g, glassColor.b, alpha)
    }
}
