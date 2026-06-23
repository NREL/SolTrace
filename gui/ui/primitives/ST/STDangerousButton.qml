import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Material
import QtQuick.Effects
import QtQuick.Layouts

import SolTrace

STButton {
    Material.foreground: Material.Red
    idle_color: App.theme.destructiveGlassColor
    down_color: Material.color(Material.Red)
}
