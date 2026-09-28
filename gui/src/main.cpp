#include "app/logging.h"

#include <QApplication>
#include <QFile>
#include <QFontDatabase>
#include <QPointer>
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
    qputenv("QT_QUICK_CONTROLS_STYLE", "Material");
    qputenv("QT_QUICK_CONTROLS_MATERIAL_VARIANT", "Dense");
    qputenv("QML_XHR_ALLOW_FILE_READ", "1");

    QApplication app(argc, argv);
    app.setOrganizationName("NLR");
    app.setOrganizationDomain("nlr.gov");
    app.setApplicationName("SolTrace");

    SolTrace::GUI::App::initialize_logging_handler();

    // Start up splash screen

    QPointer<QQmlApplicationEngine> splash_engine = new QQmlApplicationEngine();

    // If the splash fails, we exit.
    QObject::connect(
        splash_engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []() { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);

    // Load splash content
    splash_engine->load("qrc:/qt/qml/SolTrace/ui/app/Splash.qml");

    // A few spins to load up content
    app.processEvents();
    app.processEvents();


    // Initialize fonts, etc.

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

    // Main engine
    QQmlApplicationEngine engine;

    QObject::connect(
        &engine,
        &QQmlApplicationEngine::objectCreationFailed,
        &app,
        []() { QCoreApplication::exit(-1); },
        Qt::QueuedConnection);

    QObject::connect(&engine,
                     &QQmlApplicationEngine::objectCreated,
                     splash_engine,
                     &QObject::deleteLater,
                     Qt::QueuedConnection);

    // engine.

    engine.loadFromModule("SolTrace", "Main");

    return app.exec();
}
