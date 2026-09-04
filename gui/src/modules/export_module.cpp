#include "export_module.h"

#include "database/database.h"

#include <entt/entity/entity.hpp>

#include <QDir>
#include <QFile>
#include <QRandomGenerator>
#include <QRegularExpression>
#include <QSettings>
#include <QStandardPaths>
#include <QTextStream>

#include <algorithm>
#include <exception>
#include <limits>
#include <numeric>

namespace SolTrace::GUI::App {

namespace {

QString event_type_name(db::RayEventType event) {
    switch (event) {
    case db::RayEventType::CREATE: return QStringLiteral("CREATE");
    case db::RayEventType::ABSORB: return QStringLiteral("ABSORB");
    case db::RayEventType::REFLECT: return QStringLiteral("REFLECT");
    case db::RayEventType::TRANSMIT: return QStringLiteral("TRANSMIT");
    case db::RayEventType::VIRTUAL: return QStringLiteral("VIRTUAL");
    case db::RayEventType::EXIT: return QStringLiteral("EXIT");
    case db::RayEventType::UNKNOWN: return QStringLiteral("UNKNOWN");
    }

    return QStringLiteral("UNKNOWN");
}

QString csv_escape(QString value) {
    if (!value.contains(',') && !value.contains('"') && !value.contains('\n') &&
        !value.contains('\r')) {
        return value;
    }

    value.replace('"', "\"\"");
    return "\"" + value + "\"";
}

QString file_safe(QString value) {
    value = value.trimmed();
    if (value.isEmpty()) value = QStringLiteral("soltrace-result");

    value.replace(QRegularExpression(QStringLiteral("[^A-Za-z0-9._-]+")),
                  QStringLiteral("_"));
    value.replace(QRegularExpression(QStringLiteral("_+")),
                  QStringLiteral("_"));
    return value.trimmed();
}

bool ensure_directory(QString const& path) {
    if (path.isEmpty()) return false;

    QDir dir(path);
    if (dir.exists()) return true;
    return dir.mkpath(QStringLiteral("."));
}

QString path_from_url(QUrl const& url) {
    if (url.isLocalFile()) return url.toLocalFile();
    return url.toString();
}

size_t bounded_index(QRandomGenerator* random, size_t upper_exclusive) {
    if (upper_exclusive <= 1) return 0;
    if (upper_exclusive <= std::numeric_limits<quint32>::max()) {
        return random->bounded(static_cast<quint32>(upper_exclusive));
    }

    return static_cast<size_t>(random->generate64() % upper_exclusive);
}

bool write_flux_map_mesh(QString const&             path,
                         QString const&             object_name,
                         analysis::BakedFluxMapPtr const& map) {
    if (!map) return false;

    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return false;

    QTextStream out(&file);
    out.setRealNumberPrecision(10);

    out << "# SolTrace extended OBJ flux map mesh\n";
    out << "# fa <face_area>\n";
    out << "# fv <face_bin_ray_count>\n";
    out << "# power_per_ray " << map->stats.power_per_ray << '\n';
    out << "# plotted_power " << map->stats.plotted_power << '\n';
    out << "o " << object_name << '\n';

    for (auto const& vertex : map->mesh.vertex) {
        out << "v " << vertex.position.x << ' ' << vertex.position.y << ' '
            << vertex.position.z << '\n';
    }

    for (auto const& vertex : map->mesh.vertex) {
        out << "vt " << vertex.uv.x << ' ' << vertex.uv.y << '\n';
    }

    for (auto const& vertex : map->mesh.vertex) {
        out << "vn " << vertex.normal.x << ' ' << vertex.normal.y << ' '
            << vertex.normal.z << '\n';
    }

    for (qsizetype i = 0; i < map->mesh.triangles.size(); ++i) {
        auto const& triangle = map->mesh.triangles[i];

        auto const a = static_cast<qulonglong>(triangle.x) + 1;
        auto const b = static_cast<qulonglong>(triangle.y) + 1;
        auto const c = static_cast<qulonglong>(triangle.z) + 1;

        out << "f " << a << '/' << a << '/' << a << ' ' << b << '/' << b
            << '/' << b << ' ' << c << '/' << c << '/' << c << '\n';

        auto const face_area =
            i < map->face_area.size() ? map->face_area[i] : 0.0f;
        auto const face_ray_count =
            i < map->face_ray_count.size() ? map->face_ray_count[i] : 0;

        out << "fa " << face_area << '\n';
        out << "fv " << face_ray_count << '\n';
    }

    return out.status() == QTextStream::Ok;
}

} // namespace

ExportModule::ExportModule(QObject* parent) : QObject(parent) {
    QSettings settings;
    settings.beginGroup(QStringLiteral("AnalysisExport"));

    auto default_dir =
        QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);

    set_export_directory(QUrl::fromLocalFile(
        settings.value(QStringLiteral("directory"), default_dir).toString()));
    set_export_flux_map_images(
        settings.value(QStringLiteral("flux_map_images"), true).toBool());
    set_export_flux_map_meshes(
        settings.value(QStringLiteral("flux_map_meshes"), false).toBool());
    set_export_rays(settings.value(QStringLiteral("rays"), true).toBool());
    set_export_scene_copy(
        settings.value(QStringLiteral("scene_copy"), false).toBool());
    set_random_sample_rays(
        settings.value(QStringLiteral("random_sample_rays"), false).toBool());
    set_random_sample_ray_count(
        settings.value(QStringLiteral("random_sample_ray_count"), 10000)
            .toInt());
    settings.endGroup();

    connect(this, &ExportModule::export_directory_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("directory"),
                          path_from_url(export_directory()));
        update_can_export();
    });
    connect(this, &ExportModule::export_flux_map_images_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("flux_map_images"),
                          export_flux_map_images());
        update_can_export();
    });
    connect(this, &ExportModule::export_flux_map_meshes_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("flux_map_meshes"),
                          export_flux_map_meshes());
        update_can_export();
    });
    connect(this, &ExportModule::export_rays_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("rays"), export_rays());
        update_can_export();
    });
    connect(this, &ExportModule::export_scene_copy_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("scene_copy"), export_scene_copy());
        update_can_export();
    });
    connect(this, &ExportModule::random_sample_rays_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("random_sample_rays"),
                          random_sample_rays());
    });
    connect(this, &ExportModule::random_sample_ray_count_changed, this, [this] {
        QSettings settings;
        settings.beginGroup(QStringLiteral("AnalysisExport"));
        settings.setValue(QStringLiteral("random_sample_ray_count"),
                          random_sample_ray_count());
    });
}

void ExportModule::set_results(db::SimulationResultPtr results) {
    m_results = std::move(results);
    m_flux_maps.clear();

    set_current_result_name(m_results && m_results->database
                                ? m_results->database->name()
                                : QString());
    set_generated_flux_map_count(0);
    update_can_export();
}

void ExportModule::cache_flux_map(db::Entity                entity,
                                  analysis::BakedFluxMapPtr image,
                                  db::Database const*) {
    if (!m_results || !image || !entity.is_valid()) return;

    m_flux_maps.insert_or_assign(entity, std::move(image));
    set_generated_flux_map_count(static_cast<int>(m_flux_maps.size()));
}

QString ExportModule::current_result_file_stem() const {
    return file_safe(current_result_name());
}

void ExportModule::update_can_export() {
    set_can_export(m_results != nullptr &&
                   !path_from_url(export_directory()).isEmpty() &&
                   (export_flux_map_images() || export_flux_map_meshes() ||
                    export_rays() || export_scene_copy()));
}

void ExportModule::export_current() {
    if (!m_results) {
        emit notify(ANotification::warning(
            QStringLiteral("Run a simulation before exporting results.")));
        return;
    }

    auto const directory_path = path_from_url(export_directory());
    if (!ensure_directory(directory_path)) {
        emit notify(ANotification::error(QStringLiteral(
            "Could not create the export folder. Check the folder path and "
            "permissions.")));
        return;
    }

    emit notify(ANotification::info(QStringLiteral("Exporting results...")));

    QDir const directory(directory_path);
    auto const stem          = current_result_file_stem();
    int        files_written = 0;

    if (export_scene_copy()) {
        if (!m_results->database) {
            emit notify(ANotification::error(QStringLiteral(
                "Could not export the scene copy because the simulation result "
                "does not include a scene snapshot.")));
            return;
        }

        auto* database = const_cast<db::Database*>(m_results->database.get());
        auto  result   = database->export_to_simdata();

        if (!result) {
            emit notify(ANotification::error(
                QStringLiteral("Could not export the scene copy: %1")
                    .arg(result.get_failure())));
            return;
        }

        auto const scene_path =
            directory.filePath(stem + QStringLiteral("_scene.json"));

        try {
            result.get_success()->data->export_json_file(
                scene_path.toStdString());
        } catch (std::exception const& ex) {
            emit notify(ANotification::error(
                QStringLiteral("Could not write the scene copy: %1")
                    .arg(ex.what())));
            return;
        }

        ++files_written;
    }

    if (export_flux_map_images() || export_flux_map_meshes()) {
        for (auto iter = m_flux_maps.begin(); iter != m_flux_maps.end();
             ++iter) {
            auto const entity_id = QString::number(
                entt::to_integral(iter->first.value));
            auto const entity_name =
                m_results->database
                    ? file_safe(m_results->database->name_of(iter->first))
                    : QStringLiteral("entity");

            auto const base_name = stem + "_" + entity_name + "_" + entity_id;
            auto const bin_path =
                directory.filePath(base_name + QStringLiteral("_flux-map.png"));
            auto const point_path = directory.filePath(
                base_name + QStringLiteral("_ray-points.png"));
            auto const mesh_path =
                directory.filePath(base_name + QStringLiteral("_flux-mesh.obj"));

            if (export_flux_map_images()) {
                if (!iter->second->bin_map.isNull() &&
                    iter->second->bin_map.save(bin_path, "PNG")) {
                    ++files_written;
                }

                if (!iter->second->point_map.isNull() &&
                    iter->second->point_map.save(point_path, "PNG")) {
                    ++files_written;
                }
            }

            if (export_flux_map_meshes() &&
                write_flux_map_mesh(mesh_path, entity_name, iter->second)) {
                ++files_written;
            }
        }
    }

    // TODO: Thread this

    if (export_rays()) {
        auto indices = std::vector<size_t>(m_results->records.size());
        std::iota(indices.begin(), indices.end(), size_t { 0 });

        if (random_sample_rays()) {
            auto*      random = QRandomGenerator::global();
            auto const max_sample_count =
                std::min(indices.size(),
                         static_cast<size_t>(std::numeric_limits<int>::max()));
            auto const sample_count =
                std::clamp(random_sample_ray_count(),
                           0,
                           static_cast<int>(max_sample_count));
            for (size_t i = 0; i < static_cast<size_t>(sample_count); ++i) {
                auto const j = i + bounded_index(random, indices.size() - i);
                std::swap(indices[i], indices[j]);
            }
            indices.resize(static_cast<size_t>(sample_count));
            std::sort(indices.begin(), indices.end());
        }

        auto file =
            QFile(directory.filePath(stem + QStringLiteral("_rays.csv")));
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            emit notify(ANotification::error(QStringLiteral(
                "Could not write the ray data file. Check folder permissions "
                "and available disk space.")));
            return;
        }

        QTextStream out(&file);
        out << "ray_id,event_index,event_type,entity_id,entity_name,"
               "location_x,location_y,location_z,"
               "direction_x,direction_y,direction_z\n";

        for (auto const index : indices) {
            auto const& record = m_results->records[index];
            for (size_t event_index = 0; event_index < record.events.size();
                 ++event_index) {
                auto const& event = record.events[event_index];
                auto const  entity_id =
                    event.entity == entt::null
                         ? QString()
                         : QString::number(entt::to_integral(event.entity));
                auto const entity_name =
                    event.entity != entt::null && m_results->database
                        ? m_results->database->name_of(event.entity)
                        : QString();

                out << static_cast<qulonglong>(record.id) << ','
                    << static_cast<qulonglong>(event_index) << ','
                    << event_type_name(event.event) << ',' << entity_id << ','
                    << csv_escape(entity_name) << ',' << event.location.x << ','
                    << event.location.y << ',' << event.location.z << ','
                    << event.direction.x << ',' << event.direction.y << ','
                    << event.direction.z << '\n';
            }
        }

        ++files_written;
    }

    if (files_written == 0) {
        emit notify(ANotification::warning(QStringLiteral(
            "No files were exported. Enable at least one export option and "
            "try again.")));
        return;
    }

    emit notify(ANotification::info(
        QStringLiteral("Export complete: wrote %1 file(s).")
            .arg(files_written)));
}

} // namespace SolTrace::GUI::App
