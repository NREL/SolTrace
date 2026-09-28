import QtQuick
import QtCore
import QtQuick.Controls.Material

QtObject {
    id: theme

    readonly property int glassTintDefault: 0
    readonly property int glassTintRed: 1
    readonly property int glassTintPink: 2
    readonly property int glassTintPurple: 3
    readonly property int glassTintDeepPurple: 4
    readonly property int glassTintIndigo: 5
    readonly property int glassTintBlue: 6
    readonly property int glassTintLightBlue: 7
    readonly property int glassTintCyan: 8
    readonly property int glassTintTeal: 9
    readonly property int glassTintGreen: 10
    readonly property int glassTintLightGreen: 11
    readonly property int glassTintLime: 12
    readonly property int glassTintYellow: 13
    readonly property int glassTintAmber: 14
    readonly property int glassTintOrange: 15
    readonly property int glassTintDeepOrange: 16
    readonly property int glassTintBrown: 17
    readonly property int glassTintGrey: 18
    readonly property int glassTintBlueGrey: 19

    readonly property int glassClarityClear: 0
    readonly property int glassClarityBalanced: 1
    readonly property int glassClaritySolid: 2

    readonly property int fontSizeSmall: 0
    readonly property int fontSizeDefault: 1
    readonly property int fontSizeLarge: 2

    property int glassTint: glassTintDefault
    property int glassClarity: glassClarityBalanced
    property int fontSize: fontSizeDefault

    readonly property color defaultFontColor: "white"
    readonly property color fontColor: defaultFontColor
    readonly property color glassTintColor: tintColor(glassTint)
    readonly property real glassAlpha: clarityAlpha(glassClarity)
    readonly property color defaultGlassColor: Qt.rgba(0, 0, 0, 0.14)
    readonly property color glassColor: Qt.rgba(
        glassTintColor.r,
        glassTintColor.g,
        glassTintColor.b,
        glassAlpha)
    readonly property color destructiveGlassColor: Qt.rgba(
        Math.min(glassTintColor.r + 0.1, 1),
        glassTintColor.g,
        glassTintColor.b,
        glassAlpha)
    readonly property color shadedGlassColor: theme.glassColorA(
        Math.min(glassAlpha + 0.05, 1))
    readonly property color dividerColor: Qt.rgba(1, 1, 1, 0.2)

    readonly property int defaultHeaderSize: 17
    readonly property int defaultSubHeaderSize: 16
    readonly property int defaultLabelSize: 13
    readonly property int defaultComboBarTextSize: 13
    readonly property int defaultPropertyPanelHeaderSize: 14
    readonly property int defaultNormalSize: 15

    readonly property int headerSize: calculateSize(defaultHeaderSize)
    readonly property int subHeaderSize: calculateSize(defaultSubHeaderSize)
    readonly property int labelSize: calculateSize(defaultLabelSize)
    readonly property int comboBarTextSize: calculateSize(defaultComboBarTextSize)
    readonly property int propertyPanelHeaderSize: calculateSize(defaultPropertyPanelHeaderSize)
    readonly property int normalSize: calculateSize(defaultNormalSize)

    property color sunShapeGraphLineColor: "#E0D080"
    property int sunShapeGraphLineWidth: 3

    readonly property Settings settings: Settings {
        property alias glassTint: theme.glassTint
        property alias glassClarity: theme.glassClarity
        property alias fontSize: theme.fontSize
        property alias sunShapeGraphLineColor: theme.sunShapeGraphLineColor
        property alias sunShapeGraphLineWidth: theme.sunShapeGraphLineWidth
    }

    function glassColorA(alpha) {
        return Qt.rgba(glassTintColor.r, glassTintColor.g, glassTintColor.b, alpha)
    }

    function fontOffset(mode) {
        switch (mode) {
        case fontSizeSmall: return -1
        case fontSizeLarge: return 2
        default: return 0
        }
    }

    function calculateSize(baseSize) {
        return Math.max(1, Math.min(48, baseSize + fontOffset(fontSize)))
    }

    function clarityAlpha(clarity) {
        switch (clarity) {
        case glassClarityClear: return 0.01
        case glassClaritySolid: return 0.75
        default: return 0.14
        }
    }

    function tintColor(tint) {
        switch (tint) {
        case glassTintRed: return Material.color(Material.Red)
        case glassTintPink: return Material.color(Material.Pink)
        case glassTintPurple: return Material.color(Material.Purple)
        case glassTintDeepPurple: return Material.color(Material.DeepPurple)
        case glassTintIndigo: return Material.color(Material.Indigo)
        case glassTintBlue: return Material.color(Material.Blue)
        case glassTintLightBlue: return Material.color(Material.LightBlue)
        case glassTintCyan: return Material.color(Material.Cyan)
        case glassTintTeal: return Material.color(Material.Teal)
        case glassTintGreen: return Material.color(Material.Green)
        case glassTintLightGreen: return Material.color(Material.LightGreen)
        case glassTintLime: return Material.color(Material.Lime)
        case glassTintYellow: return Material.color(Material.Yellow)
        case glassTintAmber: return Material.color(Material.Amber)
        case glassTintOrange: return Material.color(Material.Orange)
        case glassTintDeepOrange: return Material.color(Material.DeepOrange)
        case glassTintBrown: return Material.color(Material.Brown)
        case glassTintGrey: return Material.color(Material.Grey)
        case glassTintBlueGrey: return Material.color(Material.BlueGrey)
        default: return Qt.rgba(0, 0, 0, 1)
        }
    }
}
