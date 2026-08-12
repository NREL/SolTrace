#pragma once

#include <QAbstractTableModel>
#include <QDebug>
#include <QPointer>

#include <span>


template <class T>
struct SMRecordMeta {
    using SMGetter = QVariant (*)(T const&);
    using SMSetter = void     (*)(T&, QVariant const&);


    const char* name   = "";
    size_t      offset = 0;
    SMGetter    getter = nullptr;
    SMSetter    setter = nullptr;
};

template <class>
struct is_shared_ptr_t : std::false_type { };
template <class U>
struct is_shared_ptr_t<std::shared_ptr<U>> : std::true_type { };
template <class T>
inline constexpr bool is_shared_ptr_v =
    is_shared_ptr_t<std::remove_cvref_t<T>>::value;

template <class T>
QVariant _sm_to_variant(T const& t) {
    return QVariant::fromValue(t);
}

template <class T>
QVariant _sm_to_variant(std::shared_ptr<T> const& t) {
    return QVariant::fromValue(t.get());
}

template <class T>
QVariant _sm_to_variant(std::unique_ptr<T> const& t) {
    return QVariant::fromValue(t.get());
}

template <class T>
QVariant _sm_to_variant(QPointer<T> const& t) {
    return QVariant::fromValue(t.data());
}

#define SM_EXPOSE_RW(MEM)                                                      \
    SMRecordMeta<Record> {                                                     \
        .name = #MEM, .offset = offsetof(Record, MEM),                         \
        .getter = [](Record const& t) -> QVariant {                            \
            return _sm_to_variant(t.MEM);                                      \
        },                                                                     \
        .setter = [](Record& t, QVariant const& a) {                           \
            using LT = std::remove_cvref_t<decltype(t.MEM)>;                   \
            t.MEM    = a.value<LT>();                                          \
        },                                                                     \
    }

#define SM_EXPOSE_RO(MEM)                                                      \
    SMRecordMeta<Record> {                                                     \
        .name = #MEM, .offset = offsetof(Record, MEM),                         \
        .getter = [](Record const& t) -> QVariant {                            \
            return _sm_to_variant(t.MEM);                                      \
        },                                                                     \
        .setter = nullptr,                                                     \
    }

// #define RECORD_META(RT, ...) \
//     template <> \
//     struct SMMetaGetter<RT> { \
//         using Record                            = RT; \
//         static inline constexpr std::array meta = { __VA_ARGS__ }; \
//     };


#define RECORD_META(RT, ...)                                                   \
    inline static constexpr auto sm_meta_getter() {                            \
        using Record              = RT;                                        \
        constexpr std::array meta = { __VA_ARGS__ };                           \
        return meta;                                                           \
    }


template <class Record>
QStringList get_header() {
    QStringList ret;
    constexpr auto meta = Record::sm_meta_getter();
    for (auto const& m : meta) {
        ret << m.name;
    }
    return ret;
}

template <class Record>
QHash<int, QByteArray> const& get_name_map() {
    constexpr auto meta = Record::sm_meta_getter();

    static QHash<int, QByteArray> ret = [&]() {
        QHash<int, QByteArray> build;

        for (int i = 0; i < std::size(meta); i++) {
            build[Qt::UserRole + i] = meta[i].name;
        }

        return build;
    }();
    return ret;
}

template <class Record>
constexpr int role_for_member_offset(size_t off) {
    constexpr auto meta = Record::sm_meta_getter();

    for (int i = 0; i < std::size(meta); i++) {
        if (meta[i].offset == off) { return Qt::UserRole + i; }
    }

    return -1;
}

#define ROLE_FOR_MEMBER(RT, MEM)                                               \
    [] {                                                                       \
        constexpr auto role = role_for_member_offset<RT>(offsetof(RT, MEM));   \
        static_assert(role != -1, "Member is not exposed as a model role");    \
        return role;                                                           \
    }()


template <class T>
constexpr bool is_qobject = std::is_base_of_v<QObject, std::remove_cvref_t<T>>;

template <class T>
struct is_shared_qobject {
    static constexpr bool value = false;
};

template <class T>
struct is_shared_qobject<std::shared_ptr<T>> {
    static constexpr bool value = is_qobject<T>;
};

template <class Record>
QVariant _record_runtime_get(Record const& r, int i) {
    constexpr auto        meta   = Record::sm_meta_getter();
    auto                  getter = meta.at(i).getter;
    if (getter) { return std::invoke(getter, r); }
    return {};
}

template <class Record>
bool _record_runtime_set(Record& r, int i, QVariant const& v) {
    constexpr auto        meta   = Record::sm_meta_getter();
    auto                  setter = meta.at(i).setter;

    if (setter) {
        std::invoke(setter, r, v);
        return true;
    }

    return false;
}

class StructTableModelBase : public QAbstractTableModel {
    Q_OBJECT
public:
    using QAbstractTableModel::QAbstractTableModel;
};

template <class Record>
class StructTableModel : public StructTableModelBase {
protected:
    QVector<Record> m_records;

    QStringList const m_header;

public:
    explicit StructTableModel(QObject* parent = nullptr)
        : StructTableModelBase(parent), m_header(get_header<Record>()) {
        static_assert(
            std::is_standard_layout_v<Record>,
            "StructTableModel Record must be standard-layout (offsetof used)");
    }

    // Header:
    QVariant headerData(int             section,
                        Qt::Orientation orientation,
                        int             role = Qt::DisplayRole) const override {
        if (orientation != Qt::Orientation::Horizontal) return {};
        if (role != Qt::DisplayRole) return {};

        return m_header.value(section);
    }

    int rowCount(QModelIndex const& parent = QModelIndex()) const override {
        if (parent.isValid()) return 0;
        return m_records.size();
    }

    int columnCount(QModelIndex const& parent = QModelIndex()) const override {
        if (parent.isValid()) return 0;
        return m_header.size();
    }

    QVariant data(QModelIndex const& index,
                  int                role = Qt::DisplayRole) const override {

        if (!index.isValid()) return {};
        if (index.row() >= m_records.size()) return {};

        auto const& item = m_records[index.row()];

        if (role == Qt::DisplayRole or role == Qt::EditRole) {
            return _record_runtime_get(item, index.column());
        }

        if (role >= Qt::UserRole) {
            auto local_role = role - Qt::UserRole;

            assert(local_role >= 0);

            if (local_role >= m_header.size()) return {};

            return _record_runtime_get(item, local_role);
        }

        return {};
    }


    bool setData(QModelIndex const& index,
                 QVariant const&    value,
                 int                role = Qt::EditRole) override {

        if (data(index, role) == value) return false;

        auto& item = m_records[index.row()];

        int location = -1;

        if (role >= Qt::UserRole) {
            location = role - Qt::UserRole;
        } else {
            location = index.column();
        }

        if (location >= m_header.size()) return false;

        bool ok = _record_runtime_set(item, location, value);

        if (!ok) return false;

        Q_EMIT dataChanged(
            index, index, { Qt::DisplayRole, Qt::EditRole, role });
        return true;
    }

    Qt::ItemFlags flags(QModelIndex const& index) const override {
        if (!index.isValid()) return Qt::NoItemFlags;

        auto const& meta = Record::sm_meta_getter();

        bool can_edit = !!meta.at(index.column()).setter;

        if (!can_edit) return Qt::ItemIsEnabled | Qt::ItemIsSelectable;

        return Qt::ItemIsEditable | Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    }

    QHash<int, QByteArray> roleNames() const override {
        static const auto roles = get_name_map<Record>();

        return roles;
    }

    bool insertRows(int                row,
                    int                count,
                    QModelIndex const& p = QModelIndex()) override {
        if (row < 0 or count <= 0) return false;

        beginInsertRows(p, row, row + count - 1);
        m_records.insert(row, count, Record {});
        endInsertRows();
        return true;
    }

    bool removeRows(int                row,
                    int                count,
                    QModelIndex const& p = QModelIndex()) override {
        if (row < 0 or count <= 0) return false;
        if (row + count > m_records.size()) return false;

        beginRemoveRows(p, row, row + count - 1);
        m_records.remove(row, count);
        endRemoveRows();
        return true;
    }

    void reset(QVector<Record> new_records = {}) {
        beginResetModel();
        m_records = new_records;
        endResetModel();
    }

    // this emits a remove signal, instead of a reset
    void remove_all() {
        if (m_records.isEmpty()) return;
        beginRemoveRows({}, 0, rowCount() - 1);
        m_records.clear();
        endRemoveRows();
    }

    Record const* get_at(int i) const {
        if (i < 0) return nullptr;
        if (i >= m_records.size()) return nullptr;
        return &m_records[i];
    }

    void append(Record const& r) {
        int rc = rowCount();
        beginInsertRows({}, rc, rc);
        m_records << r;
        endInsertRows();
    }

    void append(QVector<Record> r) {
        if (r.isEmpty()) return;
        int rc = rowCount();
        beginInsertRows({}, rc, rc + r.size() - 1);
        m_records << r;
        endInsertRows();
    }

    void replace(QVector<Record> r = {}) {
        remove_all();
        append(r);
    }

    void update(int i, Record const& r) {
        if (i < 0) return;
        if (i >= m_records.size()) return;

        m_records[i] = r;

        auto left  = index(i, 0);
        auto right = index(i, columnCount() - 1);

        Q_EMIT dataChanged(left, right);
    }

    void remove_at(int index, int count = 1) {
        if (index < 0) return;
        if (index + count > m_records.size()) return;

        beginRemoveRows(QModelIndex(), index, index + count - 1);
        m_records.remove(index, count);
        endRemoveRows();
    }

    void insert_at(int index, std::span<Record> records) {
        if (records.empty()) return;
        beginInsertRows({}, index, index + records.size() - 1);
        m_records.insert(index, records.size(), Record {});
        for (int i = 0; i < records.size(); i++) {
            m_records[index + i] = records[i];
        }
        endInsertRows();
    }


    // delete by a predicate
    template <class Function>
    void remove_by_predicate(Function&& f) {
        QVector<int> to_remove;

        for (int i = 0; i < m_records.size(); i++) {
            if (f(m_records[i])) to_remove << i;
        }

        // sort to high -> low, so indices are preserved
        std::reverse(to_remove.begin(), to_remove.end());

        for (auto i : std::as_const(to_remove)) {
            remove_at(i);
        }
    }

    // We do NOT expose a mutable view of the container
    auto const& vector() const { return m_records; }

    auto begin() const { return m_records.begin(); }
    auto end() const { return m_records.end(); }

    auto cbegin() const { return m_records.begin(); }
    auto cend() const { return m_records.end(); }
};


template <class Record>
struct StructModelAdapter : public StructTableModelBase {
    QVector<Record> m_records;

    QStringList const m_header;

protected:
    virtual QVector<Record> request_insert_blank(int row, int count) {
        return {};
    }
    virtual bool request_reset(std::span<Record>) { return false; }
    virtual bool request_remove(int i, int count) { return false; }
    virtual bool request_update(int i, Record const&) { return false; }
    virtual bool request_append(std::span<Record>) { return false; }


    // these must be infallible
    void store_push_insert(int at, std::span<Record> list) {
        this->beginInsertRows({}, at, at + list.size() - 1);
        this->m_records.insert(at, list);
        this->endInsertRows();
    }

    void store_push_append(Record r) {
        auto at = rowCount();
        this->beginInsertRows({}, at, at);
        this->m_records.insert(at, r);
        this->endInsertRows();
    }

    void store_push_remove(int at, int count) {
        this->beginRemoveRows({}, at, at + count - 1);
        this->m_records.remove(at, count);
        this->endRemoveRows();
    }
    void store_push_update(int i, Record const& item) {
        this->m_records[i] = item;

        auto col_count = this->columnCount();
        QVector<int> roles { Qt::DisplayRole, Qt::EditRole };

        for (int role = 0; role < col_count; role++) {
            roles.push_back(Qt::UserRole + role);
        }

        Q_EMIT this->dataChanged(this->index(i, 0),
                                 this->index(i, col_count - 1),
                                 roles);
    }
    void store_reset(QVector<Record> new_records = {}) {
        this->beginResetModel();
        this->m_records = new_records;
        this->endResetModel();
    }

    // this emits a remove signal, instead of a reset
    void store_remove_all() {
        if (this->m_records.empty()) { return; }

        beginRemoveRows({}, 0, this->rowCount() - 1);
        this->m_records.clear();
        this->endRemoveRows();
    }

    template <class Function>
    void store_remove_by_predicate(Function&& f) {
        QVector<int> to_remove;

        for (int i = 0; i < m_records.size(); i++) {
            if (f(m_records[i])) to_remove << i;
        }

        // sort to high -> low, so indices are preserved
        std::reverse(to_remove.begin(), to_remove.end());

        for (auto i : std::as_const(to_remove)) {
            store_push_remove(i, 1);
        }
    }

public:
    explicit StructModelAdapter(QObject* parent = nullptr)
        : StructTableModelBase(parent), m_header(get_header<Record>()) {
        static_assert(std::is_standard_layout_v<Record>,
                      "StructModelAdapter Record must be standard-layout "
                      "(offsetof used)");
    }

    // Header:
    QVariant headerData(int             section,
                        Qt::Orientation orientation,
                        int             role = Qt::DisplayRole) const override {
        if (orientation != Qt::Orientation::Horizontal) return {};
        if (role != Qt::DisplayRole) return {};

        return m_header.value(section);
    }

    int rowCount(QModelIndex const& parent = QModelIndex()) const override {
        if (parent.isValid()) return 0;
        return m_records.size();
    }

    int columnCount(QModelIndex const& parent = QModelIndex()) const override {
        if (parent.isValid()) return 0;
        return m_header.size();
    }

    QVariant data(QModelIndex const& index,
                  int                role = Qt::DisplayRole) const override {

        if (!index.isValid()) return {};
        if (index.row() >= m_records.size()) return {};

        auto const& item = m_records[index.row()];

        if (role == Qt::DisplayRole or role == Qt::EditRole) {
            return _record_runtime_get(item, index.column());
        }

        if (role >= Qt::UserRole) {
            auto local_role = role - Qt::UserRole;

            assert(local_role >= 0);

            if (local_role >= m_header.size()) return {};

            return _record_runtime_get(item, local_role);
        }

        return {};
    }


    bool setData(QModelIndex const& index,
                 QVariant const&    value,
                 int                role = Qt::EditRole) override {
        if (this->data(index, role) == value) return false;

        // COPY here
        auto item = this->m_records[index.row()];

        int location = -1;

        if (role >= Qt::UserRole) {
            location = role - Qt::UserRole;
        } else {
            location = index.column();
        }

        if (location >= this->m_header.size()) return false;

        bool ok = _record_runtime_set(item, location, value);

        if (!ok) return false;

        ok = this->request_update(index.row(), item);

        if (ok) {
            this->m_records[index.row()] = item;
            Q_EMIT this->dataChanged(
                index, index, { Qt::DisplayRole, Qt::EditRole, role });
        }

        return ok;
    }

    QHash<int, QByteArray> roleNames() const override {
        static const auto roles = get_name_map<Record>();

        return roles;
    }

    bool insertRows(int                row,
                    int                count,
                    QModelIndex const& p = QModelIndex()) override {
        if (row < 0 or count <= 0) return false;

        auto blanks = this->request_insert_blank(row, count);

        if (blanks.empty()) { return false; }

        // The virt can return a count that is DIFFERENT. we HAVE to conform.
        // Because this means the underlying structure has changed, and we need
        // to match.

        this->beginInsertRows(p, row, row + blanks.size() - 1);
        auto& l     = this->m_records;
        auto  new_l = l.mid(0, row);
        new_l.append(blanks);
        new_l.append(l.mid(row));
        this->m_records = new_l;
        this->endInsertRows();
        return true;
    }

    bool removeRows(int                row,
                    int                count,
                    QModelIndex const& p = QModelIndex()) override {
        if (row < 0 or count <= 0) return false;
        if (row + count > this->m_records.size()) return false;

        if (!request_remove(row, count)) { return false; }

        this->beginRemoveRows(p, row, row + count - 1);
        this->m_records.remove(row, count);
        this->endRemoveRows();
        return true;
    }


    bool append(Record const& r) {
        if (!this->request_append({ r })) { return false; }

        int rc = this->rowCount();
        this->beginInsertRows({}, rc, rc);
        this->m_records << r;
        this->endInsertRows();

        return true;
    }

    bool append(QVector<Record> r) {
        if (r.isEmpty()) return false;
        if (!this->request_append(r)) { return false; }

        int rc = this->rowCount();
        this->beginInsertRows({}, rc, rc + r.size() - 1);
        this->m_records << r;
        this->endInsertRows();
        return true;
    }

    // We do NOT expose a mutable view of the container
    auto const& vector() const { return this->m_records; }

    auto begin() const { return this->m_records.begin(); }
    auto end() const { return this->m_records.end(); }

    auto cbegin() const { return this->m_records.begin(); }
    auto cend() const { return this->m_records.end(); }

    Record const* get_at(int i) const {
        if (i < 0) return nullptr;
        if (i >= m_records.size()) return nullptr;
        return &m_records[i];
    }
};
