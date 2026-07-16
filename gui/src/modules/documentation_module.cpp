#include "documentation_module.h"
#include "utilities/qt_file_helpers.h"
#include <QDebug>
#include <QDir>
#include <QFileInfo>
#include <QFileInfoList>
#include <QStringList>

#define DOCS_LOAD_INFO false
#define DOCS_DEBUG_INFO false

namespace SolTrace::GUI::App {

DocumentationModule::DocumentationModule(QObject* parent)
    : QObject(parent),
      m_directory_path(":/docs"),
      m_status(new StatusComponent(this)) {
    connect(this,
            &DocumentationModule::locale_changed,
            this,
            &DocumentationModule::load);

    load();
}

void DocumentationModule::load() {
    if (DOCS_LOAD_INFO) {
        qDebug() << "eng map size: " << m_docs[Locale::EN].size();
        qDebug() << "esp map size: " << m_docs[Locale::ES].size();
    }

    if (m_locale == Locale::EN && m_docs.contains(Locale::EN) &&
        !m_docs[Locale::EN].empty())
        return;
    if (m_locale == Locale::ES && m_docs.contains(Locale::ES) &&
        !m_docs[Locale::ES].empty())
        return;

    if (DOCS_LOAD_INFO) qDebug() << "Loading docs...";

    doc_walker(locale_directory());

    if (DOCS_LOAD_INFO) {
        for (auto locale = m_docs.begin(); locale != m_docs.end(); locale++) {
            QString name = (locale.key() == Locale::EN) ? "English" : "Spanish";
            qDebug() << "==========[" + name + " DOCS]==========";
            for (auto doc = locale.value().begin(); doc != locale.value().end();
                 doc++) {
                qDebug() << doc.key() << " " << doc.value()->metadata("title")
                         << " ";
            }
        }

        qDebug() << "Finished loading docs...";
    }

    set_version(version() + 1);
}

void DocumentationModule::doc_walker(const QString& dir_path,
                                     const QString& key_prefix,
                                     int            depth) {
    QDir dir(dir_path);

    if (!dir.exists()) {
        if (DOCS_DEBUG_INFO)
            qDebug() << "Directory does not exist: " << dir_path;
        return;
    }

    QString dir_name = dir.dirName();

    // Process files in this directory
    QFileInfoList file_infos = dir.entryInfoList({ "*.md" }, QDir::Files);

    for (const QFileInfo& file_info : file_infos) {
        QFile file(file_info.canonicalFilePath());

        QString name = file_info.fileName().replace(".md", "");

        QString new_key;

        if (depth == 0) {
            new_key = name;
        } else if (name == dir_name) {
            new_key = key_prefix;
        } else {
            new_key = key_prefix.isEmpty() ? name : key_prefix + "." + name;
        }

        m_docs[m_locale].insert(new_key, parse_markdown_file(file));
    }

    QFileInfoList subdirs =
        dir.entryInfoList(QDir::Dirs | QDir::NoDotAndDotDot);

    for (const QFileInfo& subdir_info : subdirs) {
        QString name = subdir_info.fileName();
        QString new_prefix =
            key_prefix.isEmpty() ? name : key_prefix + "." + name;

        doc_walker(subdir_info.absoluteFilePath(), new_prefix, depth + 1);
    }
}

QString DocumentationModule::get(QString key, QString metadata_key) {
    auto* result = m_docs[m_locale].value(key, nullptr);

    if (result == nullptr) return "Error: invalid doc key " + key;

    if (metadata_key.isEmpty()) return result->body();

    return result->metadata(metadata_key);
}

QString DocumentationModule::locale_string() {
    switch (m_locale) {
    case Locale::EN: return "en";
    case Locale::ES: return "es";
    default: return "";
    }
}

QString DocumentationModule::locale_name_string() {
    switch (m_locale) {
    case Locale::EN: return "English";
    case Locale::ES: return "Spanish";
    default: return "";
    }
}

QString DocumentationModule::locale_directory() {
    return m_directory_path + "/" + locale_string();
}

} // namespace SolTrace::GUI::App
