#include "modules/database_module.h"

#include "simulation_data.hpp"

#include <gtest/gtest.h>

#include <QFileInfo>
#include <QTemporaryDir>
#include <QUrl>

#include <filesystem>
#include <vector>

namespace {

namespace SD = SolTrace::Data;
namespace App = SolTrace::GUI::App;

std::filesystem::path simple_test_case_path() {
    return std::filesystem::path(SOLTRACE_REPO_ROOT) /
           "gui/assets/examples/simple_test_case.stinput";
}

} // namespace

TEST(DatabaseModuleSave, SaveCurrentWritesLoadableSceneFileAndNotifies) {
    App::DatabaseModule module;
    std::vector<ANotification> notifications;

    QObject::connect(&module,
                     &App::DatabaseModule::notify,
                     [&notifications](ANotification notification) {
                         notifications.push_back(notification);
                     });

    module.load_url(QUrl(), QStringLiteral("save-test"));

    ASSERT_EQ(module.rowCount(), 1);
    ASSERT_NE(module.current_database(), nullptr);

    SD::SimulationData source;
    ASSERT_TRUE(source.import_from_file(simple_test_case_path().string()));
    module.current_database()->import(source);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());

    auto const save_path = temp_dir.filePath(QStringLiteral("saved_scene.json"));
    module.save_current(QUrl::fromLocalFile(save_path));

    QFileInfo const saved_file(save_path);
    ASSERT_TRUE(saved_file.isFile());
    EXPECT_GT(saved_file.size(), 0);

    SD::SimulationData saved;
    ASSERT_NO_THROW(saved.import_json_file(save_path.toStdString()));
    EXPECT_EQ(saved.get_number_of_ray_sources(),
              source.get_number_of_ray_sources());
    EXPECT_EQ(saved.get_number_of_elements(), source.get_number_of_elements());

    ASSERT_FALSE(notifications.empty());
    auto const& last_notification = notifications.back();
    EXPECT_EQ(last_notification.type, ANotification::INFO);
    EXPECT_EQ(last_notification.message,
              QStringLiteral("File successfully saved."));
}
