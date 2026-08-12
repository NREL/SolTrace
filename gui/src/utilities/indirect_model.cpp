#include "indirect_model.h"


IndirectTableModel::IndirectTableModel(QObject* parent)
    : QAbstractTableModel(parent) { }

IndirectTableModel::~IndirectTableModel() = default;

void IndirectTableModel::add_properties(QVector<IProperty> lps) {

    for (auto const& p : lps) {
        if (std::string_view(p.display_name).empty()) {
            throw std::invalid_argument("All display names must be non-empty");
        }

        if (!p.getter) {
            throw std::invalid_argument(
                "All properties must have valid getters");
        }

        auto property_index = m_header.size();
        m_header.push_back(p.display_name);

        m_properties.push_back(p);


        m_name_map[Qt::UserRole + property_index] =
            m_properties[property_index].display_name;
    }
}

QVariant IndirectTableModel::headerData(int             section,
                                        Qt::Orientation orientation,
                                        int             role) const {
    if (orientation != Qt::Orientation::Horizontal) return {};
    if (role != Qt::DisplayRole) return {};

    return m_header.value(section);
}

int IndirectTableModel::rowCount(const QModelIndex& parent) const {
    if (parent.isValid()) return 0;
    return _record_count();
}

int IndirectTableModel::columnCount(const QModelIndex& parent) const {
    if (parent.isValid()) return 0;
    return m_header.size();
}

QVariant IndirectTableModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid()) return {};
    if (index.row() >= _record_count()) return {};

    auto resolved_index = index.row();

    if (role == Qt::DisplayRole or role == Qt::EditRole) {
        return m_properties.value(index.column()).getter(resolved_index);
    }

    if (role >= Qt::UserRole) {
        auto local_role = role - Qt::UserRole;

        if (local_role >= m_header.size() or local_role < 0) return {};

        return m_properties.value(local_role).getter(resolved_index);
    }

    return {};
}

bool IndirectTableModel::setData(const QModelIndex& index,
                                 const QVariant&    value,
                                 int                role) {
    if (data(index, role) == value) return false;

    auto resolved_index = index.row();

    int location = -1;

    if (role >= Qt::UserRole) {
        location = role - Qt::UserRole;
    } else {
        location = index.column();
    }

    if (location >= m_header.size()) return false;

    auto const& setter = m_properties.value(location).setter;

    if (!setter) return false;

    bool ok = setter(resolved_index, value);

    if (ok) {
        Q_EMIT dataChanged(
            index, index, { Qt::DisplayRole, Qt::EditRole, role });
    }

    return ok;
}

Qt::ItemFlags IndirectTableModel::flags(const QModelIndex& index) const {
    if (!index.isValid()) return Qt::NoItemFlags;

    bool can_edit = !!m_properties.value(index.column()).setter;

    if (!can_edit) return Qt::ItemIsEnabled | Qt::ItemIsSelectable;

    return Qt::ItemIsEditable | Qt::ItemIsSelectable | Qt::ItemIsEnabled;
}

QHash<int, QByteArray> IndirectTableModel::roleNames() const {
    return m_name_map;
}

bool IndirectTableModel::removeRows(int row, int count, const QModelIndex& p) {
    if (row < 0 or count <= 0) return false;
    if ((row + count) > _record_count()) return false;


    if (!_can_delete_at(row, count)) return false;

    beginRemoveRows(p, row, row + count - 1);
    _delete_at(row, count);
    endRemoveRows();
    return true;
}

bool IndirectTableModel::ask_append_record(QVariant data) {
    auto can_add_new_value = this->_can_append_new(data);

    if (!can_add_new_value) { return false; }

    auto at = _record_count();

    beginInsertRows(QModelIndex(), at, at);
    _append_new(data);
    endInsertRows();
    return true;
}

void IndirectTableModel::reset(QList<QVariant> new_records) {
    beginResetModel();
    _clear();
    for (auto const& v : std::as_const(new_records)) {
        if (_can_append_new(v)) { _append_new(v); }
    }
    endResetModel();
}

void IndirectTableModel::remove_all() {
    if (_record_count() == 0) return;
    beginRemoveRows(QModelIndex(), 0, rowCount() - 1);
    _clear();
    endRemoveRows();
}

void IndirectTableModel::notify_update(int i) {
    if (i < 0) return;
    if (i >= _record_count()) return;

    auto left  = index(i, 0);
    auto right = index(i, columnCount() - 1);

    Q_EMIT dataChanged(left, right);
}
