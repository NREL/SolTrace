#pragma once

#include <QAbstractTableModel>

struct IProperty {
    const char*                                display_name = "";
    std::function<QVariant(size_t index)>             getter;
    std::function<bool(size_t index, QVariant value)> setter;
};

// TODO: More caching

class IndirectTableModel : public QAbstractTableModel {
    QVector<IProperty>         m_properties;
    QVector<QString>           m_header;
    QHash<int, QByteArray>     m_name_map;

protected:
    virtual bool _can_append_new(QVariant const&) { return true; }
    virtual void _append_new(QVariant)            = 0;
    virtual bool _can_delete_at(size_t, size_t) { return true; }
    virtual void _delete_at(size_t, size_t)       = 0;
    virtual int  _record_count() const            = 0;
    virtual void _clear()                         = 0;

    void add_properties(QVector<IProperty>);

public:
    IndirectTableModel(QObject* parent = nullptr);
    ~IndirectTableModel();

public: // Model implementations
    QVariant headerData(int             section,
                        Qt::Orientation orientation,
                        int             role = Qt::DisplayRole) const override;

    int rowCount(QModelIndex const& parent = QModelIndex()) const override;

    int columnCount(QModelIndex const& parent = QModelIndex()) const override;

    QVariant data(QModelIndex const& index,
                  int                role = Qt::DisplayRole) const override;

    bool setData(QModelIndex const& index,
                 QVariant const&    value,
                 int                role = Qt::EditRole) override;

    Qt::ItemFlags flags(QModelIndex const& index) const override;

    QHash<int, QByteArray> roleNames() const override;


    bool removeRows(int                row,
                    int                count,
                    QModelIndex const& p = QModelIndex()) override;

public slots:

    bool ask_append_record(QVariant data);

    void reset(QList<QVariant> new_records = {});

    void remove_all();

    void notify_update(int i);
};
