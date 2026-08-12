Project(GUIUnitTests)

find_package(Qt6 6.11 REQUIRED COMPONENTS
    Core Gui Quick Quick3D Widgets Graphs Concurrent
)

set(GUI_UNIT_TEST_SRC
    ../../../gui/src/database/entity.h
    ../../../gui/src/database/components.cpp
    ../../../gui/src/database/conversion.cpp
    ../../../gui/src/database/database.cpp
    ../../../gui/src/database/database_clone.cpp
    ../../../gui/src/database/database_import.cpp
    ../../../gui/src/database/database_export.cpp
    ../../../gui/src/database/database_notification.cpp
    ../../../gui/src/modules/database_module.cpp
    ../../../gui/src/utilities/asynctask.cpp
    ../../../gui/src/utilities/euler_angles.cpp
    ../../../gui/src/utilities/euler_angles.h
    ../../../gui/src/utilities/notification.h
    ../../../gui/src/utilities/structmodel.h
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
