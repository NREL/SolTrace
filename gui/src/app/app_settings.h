#pragma once

namespace SolTrace::GUI::App {

class AppData;

class AppSettings {
public:
    static void load_session(AppData& app);
    static void save_session(AppData const& app);
    static void clear_session();
};

} // namespace SolTrace::GUI::App
