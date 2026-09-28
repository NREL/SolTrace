Project(GUIUnitTests)

find_package(Qt6 6.11 REQUIRED COMPONENTS
    Core Gui Quick Quick3D Widgets Graphs Concurrent
)

set(GUI_UNIT_TEST_SRC
    ../../../gui/src/data/entity.h
    ../../../gui/src/data/components.cpp
    ../../../gui/src/data_io/conversion.cpp
    ../../../gui/src/data/database.cpp
    ../../../gui/src/data/database_clone.cpp
    ../../../gui/src/data_io/database_import.cpp
    ../../../gui/src/data_io/database_export.cpp
    ../../../gui/src/data_io/database_notification.cpp
    ../../../gui/src/modules/database_module.cpp
    ../../../gui/src/support/asynctask.cpp
    ../../../gui/src/support/euler_angles.cpp
    ../../../gui/src/support/euler_angles.h
    ../../../gui/src/support/notification.h
    ../../../gui/src/support/structmodel.h
    conversion_test.cpp
    database_module_test.cpp
    database_roundtrip_test.cpp
    euler_angles_test.cpp
)

add_executable(GUIUnitTests
    ${GUI_UNIT_TEST_SRC}
)

target_include_directories(GUIUnitTests
    PRIVATE
        ../../../gui/src
        ../../../gui/src/data
        ../../../gui/src/data_io
        ../../../gui/src/modules
        ../../../gui/src/support
        ../../../coretrace/simulation_data
        ../../../coretrace/simulation_results
        ../../../coretrace/simulation_runner
        ../../../coretrace/simulation_runner/native_runner
)

target_link_libraries(GUIUnitTests
    PRIVATE
        Qt6::Core
        Qt6::Concurrent
        Qt6::Gui
        Qt6::Qml
        Qt6::Quick
        EnTT
        magic_enum
        simdata
        simresult
        native_runner
        GTest::gtest_main
)

target_compile_definitions(GUIUnitTests
    PRIVATE
        SOLTRACE_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/../../.."
)

set_target_properties(GUIUnitTests PROPERTIES
    AUTOMOC ON
    CXX_STANDARD 20
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(GUIUnitTests PRIVATE --coverage)
    target_link_libraries(GUIUnitTests PRIVATE gcov)
endif()

if(UNIX)
    target_link_libraries(GUIUnitTests PRIVATE pthread)
endif()

add_test(NAME GUIUnitTests COMMAND GUIUnitTests)
