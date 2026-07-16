cmake_minimum_required(VERSION 3.16)

project(SolTraceApp VERSION 4.0.0 LANGUAGES CXX)

include(cmake/CPM.cmake)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Git QUIET)

set(SOLTRACE_GIT_COMMIT "unknown")
set(SOLTRACE_GIT_BRANCH "unknown")
set(SOLTRACE_GIT_TAG "")
set(SOLTRACE_GIT_DESCRIBE "unknown")
set(SOLTRACE_GIT_DIRTY "unknown")
set(SOLTRACE_IS_PRERELEASE ON)

if(Git_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE SOLTRACE_GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE SOLTRACE_GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --exact-match HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE SOLTRACE_GIT_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE SOLTRACE_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} status --porcelain
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE SOLTRACE_GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )

    if(SOLTRACE_GIT_STATUS STREQUAL "")
        set(SOLTRACE_GIT_DIRTY "false")
    else()
        set(SOLTRACE_GIT_DIRTY "true")
    endif()
endif()

if(SOLTRACE_GIT_COMMIT STREQUAL "")
    set(SOLTRACE_GIT_COMMIT "unknown")
endif()
if(SOLTRACE_GIT_BRANCH STREQUAL "")
    set(SOLTRACE_GIT_BRANCH "unknown")
endif()
if(SOLTRACE_GIT_DESCRIBE STREQUAL "")
    set(SOLTRACE_GIT_DESCRIBE "unknown")
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/app_build_info.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/app_build_info.h"
    @ONLY
)

find_package(Qt6 6.11 REQUIRED COMPONENTS
    Core Gui Qml Quick Quick3D Widgets Graphs Concurrent
)

set(QT_QML_GENERATE_QMLLS_INI ON)

qt_standard_project_setup(REQUIRES 6.11)


# Setup icons
if(APPLE)
    set(MACOSX_BUNDLE_ICON_FILE SolTrace.icns)
    set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/appicon/SolTrace.icns")
        set_source_files_properties(${app_icon_macos} PROPERTIES
               MACOSX_PACKAGE_LOCATION "Resources")

    qt_add_executable(SolTrace MACOSX_BUNDLE src/main.cpp ${app_icon_macos})
else()
    qt_add_executable(SolTrace
        src/main.cpp
    )
endif()



# Add all QML files

file(GLOB_RECURSE SOLTRACE_QML_FILES
    CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "ui/*.qml"
)
list(SORT SOLTRACE_QML_FILES)

# Add all cpp files

file(GLOB_RECURSE SOLTRACE_CPP_FILES
    CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "src/*.cpp"
    "src/*.h"
)
list(REMOVE_ITEM SOLTRACE_CPP_FILES "src/main.cpp")
list(SORT SOLTRACE_CPP_FILES)

# TODO tighten bounds on files
file(GLOB_RECURSE SOLTRACE_ASSET_FILES
    CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "assets/*"
)
list(FILTER SOLTRACE_ASSET_FILES EXCLUDE REGEX "/\\.DS_Store$")
list(SORT SOLTRACE_ASSET_FILES)

# Add all documents
file(GLOB_RECURSE SOLTRACE_DOC_FILES
    CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "docs/*.md"
)
list(SORT SOLTRACE_DOC_FILES)

set_source_files_properties(ui/App.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)

qt_add_resources(SolTrace "assets"
    PREFIX "/"
    FILES
        ${SOLTRACE_ASSET_FILES}
)

qt_add_resources(SolTrace "docs"
    PREFIX "/"
    FILES
        ${SOLTRACE_DOC_FILES}
)

qt_add_qml_module(SolTrace
    URI SolTrace
    VERSION 1.0

    QML_FILES
        ${SOLTRACE_QML_FILES}

    SOURCES
        ${SOLTRACE_CPP_FILES}

    DEPENDENCIES
        QtQuick
        QtQuick3D
)

# Add examples folder
if(APPLE)
    # Copies into .app/Contents/Resources/examples
    file(GLOB EXAMPLE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/assets/examples/*")
    list(FILTER EXAMPLE_FILES EXCLUDE REGEX "/\\.DS_Store$")
    set_source_files_properties(${EXAMPLE_FILES} PROPERTIES
        MACOSX_PACKAGE_LOCATION "Resources/examples")
    target_sources(SolTrace PRIVATE ${EXAMPLE_FILES})
else()
    # Copies alongside the binary
    add_custom_command(TARGET SolTrace POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/assets/examples"
        "$<TARGET_FILE_DIR:SolTrace>/examples")
endif()

CPMAddPackage(
    NAME EnTT
    VERSION 3.16.0
    GITHUB_REPOSITORY skypjack/entt
    DOWNLOAD_ONLY True
)

if(EnTT_ADDED)
    add_library(EnTT INTERFACE)
    target_include_directories(EnTT SYSTEM INTERFACE ${EnTT_SOURCE_DIR}/src)
endif()

CPMAddPackage(
    NAME magic_enum
    GITHUB_REPOSITORY Neargye/magic_enum
    GIT_TAG v0.9.7
)

target_link_libraries(SolTrace PRIVATE EnTT magic_enum)

if(TARGET embree_runner)
    target_compile_definitions(SolTrace PRIVATE SOLTRACE_HAS_EMBREE_RUNNER)
    target_link_libraries(SolTrace PRIVATE embree_runner)
endif()

if(TARGET optix_runner)
    target_compile_definitions(SolTrace PRIVATE SOLTRACE_HAS_OPTIX_RUNNER)
    target_link_libraries(SolTrace PRIVATE optix_runner)
endif()

target_include_directories(SolTrace PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(SolTrace PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_include_directories(SolTrace PRIVATE ../coretrace/simulation_data)
target_include_directories(SolTrace PRIVATE ../coretrace/simulation_results)
target_include_directories(SolTrace PRIVATE ../coretrace/simulation_runner)

set_target_properties(SolTrace PROPERTIES CXX_STANDARD 20)

target_compile_definitions(SolTrace PRIVATE QT_MESSAGELOGCONTEXT)

set_target_properties(SolTrace PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER "gov.nlr.SolTrace"
    MACOSX_BUNDLE_BUNDLE_NAME "SolTrace"
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

if(UNIX AND NOT APPLE)
    # Avoid colliding with the build-tree QML module directory for URI SolTrace.
    set_target_properties(SolTrace PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
    )
endif()

target_link_libraries(SolTrace
    PRIVATE
        Qt6::Core
        Qt6::Qml
        Qt6::Gui
        Qt6::Quick
        Qt6::Quick3D
        Qt6::Widgets
        Qt6::Graphs
        Qt6::Concurrent
        simdata
        simresult
        native_runner
)

include(GNUInstallDirs)
install(TARGETS SolTrace
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/scripts"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/SolTrace/assets"
    PATTERN ".DS_Store" EXCLUDE
)

qt_generate_deploy_qml_app_script(
    TARGET SolTrace
    OUTPUT_SCRIPT deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
    DEPLOY_USER_QML_MODULES_ON_UNSUPPORTED_PLATFORM
)
install(SCRIPT ${deploy_script})
