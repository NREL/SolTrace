#include "module_common.h"

namespace SolTrace::GUI::App {

StatusComponent::StatusComponent(QObject *parent) :
    QObject(parent)
{}

void StatusComponent::mark_incomplete(const QString &message)
{
    m_status = Status::Incomplete;
    m_message = message;
}

void StatusComponent::mark_loading(const QString &message)
{
    m_status = Status::Loading;
    m_message = message;
}

void StatusComponent::mark_ready(const QString &message)
{
    m_status = Status::Ready;
    m_message = message;
}

void StatusComponent::mark_error(const QString &message)
{
    m_status = Status::Error;
    m_message = message;
}

void StatusComponent::mark_stale(const QString &message)
{
    m_status = Status::Stale;
    m_message = message;
}

void StatusComponent::mark_complete(const QString &message)
{
    m_status = Status::Complete;
    m_message = message;
}

void StatusComponent::mark_unset(const QString &message)
{
    m_status = Status::Unset;
    m_message = message;
}


QHash<int, QByteArray> PresetComponentBase::roleNames() const
{
    return {
        { NameRole, "name"},
        { DescriptionRole, "description"}
    };
}

bool PresetComponentBase::create(const QString &name, const QString &description)
{
    return create_impl(name, description);
}

bool PresetComponentBase::modify(const QString& name) {
    return modify_impl(name);
}

bool PresetComponentBase::load(const QString &name)
{
    return load_impl(name);
}

bool PresetComponentBase::save(const QString &name, const QString &description)
{
    return save_impl(name, description);
}

bool PresetComponentBase::remove(const QString &name)
{
    return remove_impl(name);
}

bool PresetComponentBase::import_preset(const QString &filepath)
{
    return import_preset_impl(filepath);
}

bool PresetComponentBase::export_preset(const QString &name, const QString &filepath)
{
    return export_preset_impl(name, filepath);

}


} // namespace SolTrace::GUI::App
