#pragma once

#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileInfoList>
#include <QMap>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
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
      : m_metadata(metadata), m_body(body), m_blocks_json("[]") {}

  MarkdownDocument(const QStringList &metadata,
                   const QString &body,
                   const QString &blocks_json)
      : m_metadata(metadata), m_body(body), m_blocks_json(blocks_json) {}

  QString metadata(QString key) {
    int index = m_metadata.indexOf(key.toLower());
    if (index == -1 || m_metadata.count() <= index + 1) {
      return "Error: metadata " + key + " not found";
    }
    return m_metadata.at(index + 1);
  }

  QString body() { return m_body; }

  QString blocks_json() { return m_blocks_json; }

  void set_metadata(const QString &key, const QString &value) {
    m_metadata.append(key);
    m_metadata.append(value);
  }

private:
  QStringList m_metadata;
  const QString m_body;
  const QString m_blocks_json;
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

inline MarkdownDocument *parse_processed_doc_file(QFile &file) {
  QFileInfo info(file);
  QString raw = read_file(file);

  QJsonParseError error;
  auto document = QJsonDocument::fromJson(raw.toUtf8(), &error);
  if (error.error != QJsonParseError::NoError || !document.isObject()) {
    qDebug() << "Cannot parse processed docs:" << info.canonicalFilePath()
             << error.errorString();
    return new MarkdownDocument({}, "ERROR: invalid processed documentation");
  }

  auto object = document.object();
  auto metadata_object = object.value("metadata").toObject();
  QStringList metadata;

  for (auto iter = metadata_object.begin(); iter != metadata_object.end(); ++iter) {
    metadata.append(iter.key());
    metadata.append(iter.value().toString());
  }

  auto blocks = object.value("blocks").toArray();
  QStringList body_parts;
  for (auto const &value : blocks) {
    auto block = value.toObject();
    if (block.value("type").toString() == "text") {
      body_parts.append(block.value("content").toString());
    }
  }

  auto blocks_json = QString::fromUtf8(QJsonDocument(blocks).toJson(QJsonDocument::Compact));
  return new MarkdownDocument(metadata, body_parts.join(""), blocks_json);
};

} // namespace SolTrace::GUI::App
