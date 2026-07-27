#include "logging.h"

#include <QApplication>
#include <QFile>
#include <QFontDatabase>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QQuickWindow>

/// List of known fonts
constexpr auto font_list = std::array {
    ":/assets/fonts/computer-modern/cmunrm.ttf",
    ":/assets/fonts/computer-modern/cmunbx.ttf",
    ":/assets/fonts/computer-modern/cmunti.ttf",
    ":/assets/fonts/computer-modern/cmunbi.ttf",
    ":/assets/fonts/roboto/Roboto-Regular.ttf",
    ":/assets/fonts/roboto/Roboto-Italic.ttf",
    ":/assets/fonts/roboto/Roboto-BoldItalic.ttf",
    ":/assets/fonts/roboto/Roboto-Bold.ttf",
    ":/assets/fonts/font-awesome/fa_solid_7.otf",
};

int main(int argc, char* argv[]) {
    qputenv("QT_QUICK_CONTROLS_MATERIAL_VARIANT", "Dense");
    qputenv("QML_XHR_ALLOW_FILE_READ", "1");

    QApplication app(argc, argv);
    app.setOrganizationName("NLR");
    app.setOrganizationDomain("nlr.gov");
    app.setApplicationName("SolTrace");

    SolTrace::GUI::App::initialize_logging_handler();

    /*
    #ifdef QT_QML_DEBUG
        QQmlDebuggingEnabler::startTcpDebugServer(
            3768,
            QQmlDebuggingEnabler::DoNotWaitForClient,
            // QQmlDebuggingEnabler::WaitForClient,
            QStringLiteral("127.0.0.1"));
    #endif
    */

    // Load fonts
    for (auto font : font_list) {
        auto result = QFontDatabase::addApplicationFont(font);
        if (result < 0) { qWarning() << "Unable to load font:" << font; }
    }

    QQmlApplicationEngine engine;

    QObject::connect(
        &engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []() { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);

    engine.loadFromModule("SolTrace", "Main");

    return app.exec();
}
