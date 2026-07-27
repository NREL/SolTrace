import QtCore
import QtQuick
import QtQuick.Dialogs

import SolTrace

Item {
    id: root

    visible: false
    width: 0
    height: 0

    property alias recent_files: file_settings.recent_files

    function load_new() {
        App.file_source.load_new()
    }

    function open_file() {
        if (App.file_source.open_file_dialog()) {
            return
        }

        if (file_settings.last_selected_file.toString() !== "")
            openFileDialog.selectedFile = file_settings.last_selected_file
        openFileDialog.open()
    }

    function load_recent(file_path) {
        var file_url = file_url_text(file_path)
        add_files(file_url)
        App.file_source.load_url(Qt.url(file_url))
    }

    function file_url_text(file_path) {
        return String(file_path)
    }

    function file_name(file_path) {
        var file_url = file_url_text(file_path)
        return decodeURIComponent(file_url.split('/').pop())
    }

    function open_example() {
        exampleFileDialog.open()
    }

    function recent_files_array() {
        var files = []

        if (!file_settings.recent_files) {
            return files
        }

        for (var i = 0; i < file_settings.recent_files.length; ++i) {
            files.push(file_url_text(file_settings.recent_files[i]))
        }

        return files
    }

    function set_recent_files(files) {
        file_settings.recent_files = files
    }

    function clear_recent_files() {
        set_recent_files([])
    }

    function add_files(file_path) {
        var normalized_file_path = file_url_text(file_path)
        var files = recent_files_array()
        var index = files.indexOf(normalized_file_path)

        if (index !== -1) {
            files.splice(index, 1)
        }

        files.unshift(normalized_file_path)

        if (files.length > 5) {
            files.pop()
        }

        // Assign a fresh array so views using recent_files are notified.
        set_recent_files(files)
    }

    QtObject {
        id: file_settings

        property var recent_files: []
        property url last_selected_file
        property url last_selected_folder: StandardPaths.standardLocations(
                                               StandardPaths.DocumentsLocation)[0]
    }

    Settings {
        id: file_settings_storage

        category: "file_history"

        property alias recent_files: file_settings.recent_files
        property alias last_selected_file: file_settings.last_selected_file
        property alias last_selected_folder: file_settings.last_selected_folder
    }

    FileDialog {
        id: openFileDialog

        currentFolder: file_settings.last_selected_folder

        onAccepted: {
            var str_file = String(selectedFile)

            file_settings.last_selected_file = selectedFile
            file_settings.last_selected_folder = str_file.substring(
                        0, str_file.lastIndexOf("/"))
            root.add_files(selectedFile.toString())
            App.file_source.load_url(selectedFile)
        }
    }


    FileDialog {
        id: exampleFileDialog

        currentFolder: App.file_source.examples_folder

        onAccepted: {
            var str_file = String(selectedFile)

            file_settings.last_selected_file = selectedFile
            file_settings.last_selected_folder = str_file.substring(
                        0, str_file.lastIndexOf("/"))
            root.add_files(selectedFile.toString())
            App.file_source.load_url(selectedFile)
        }
    }
}
