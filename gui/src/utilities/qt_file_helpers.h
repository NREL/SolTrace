#pragma once

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileInfoList>
#include <QMap>
#include <QStringList>
#include <QTextStream>

namespace SolTrace::GUI::App {

QFile inline load_file(const QDir &dir, const QString &filename) {
  return QFile(dir.filePath(filename));
}

QString inline read_file(QFile &file) {
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug() << "Cannot open file: " << file.fileName();
    return {};
  }
  QTextStream in(&file);
  return in.readAll();
}

QStringList inline get_subfolders(const QDir &dir) {
  QStringList full_paths;
  const QStringList names = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
  full_paths.reserve(names.size());
  for (const QString &name : names)
    full_paths << dir.absolutePath() + "/" + name;
  return full_paths;
}

class MarkdownDocument {
public:
  MarkdownDocument(const QStringList &metadata, const QString &body)
      : m_metadata(metadata), m_body(body) {}

  QString metadata(QString key) {
    int index = m_metadata.indexOf(key.toLower());
    if (index == -1 || m_metadata.count() <= index + 1) {
      return "Error: metadata " + key + " not found";
    }
    return m_metadata.at(index + 1);
  }

  QString body() { return m_body; }

  void set_metadata(const QString &key, const QString &value) {
    m_metadata.append(key);
    m_metadata.append(value);
  }

private:
  QStringList m_metadata;
  const QString m_body;
};

inline MarkdownDocument  *parse_markdown_file(QFile &file) {
  QFileInfo info(file);
  QString raw = read_file(file);

  QStringList metadata;
  QString body;

  // Case 1: No frontmatter
  if (!raw.trimmed().startsWith("---")) {
    body = raw.trimmed();
  }

  // Case 2: Invalid frontmatter block
  auto end = raw.indexOf("---", 3);

  if (end == -1) {
    qDebug() << "Unclosed frontmatter block found in " +
                    info.canonicalFilePath();
    body = "ERROR: Unclosed frontmatter block";
  }

  // Case 3: Parse frontmatter
  auto frontmatter = raw.mid(3, end - 3).trimmed();

  for (const QString &line : frontmatter.trimmed().split("\n")) {
    auto colon = line.indexOf(":");
    QString key = line.left(colon).trimmed();
    QString value = line.mid(colon + 1).trimmed().replace("\"", "");

    metadata.append(key);
    metadata.append(value);
  }

  body = raw.mid(end + 3).trimmed();

  return new MarkdownDocument(metadata, body);
};

} // namespace SolTrace::GUI::App
