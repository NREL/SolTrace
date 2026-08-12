#include "database_clone.h"

#include "components.h"

namespace db {

/// Remap entities on clone
struct EntityMapper {
    entt::registry&                                new_registry;
    std::unordered_map<entt::entity, entt::entity> old_to_new_map;

    entt::entity operator()(entt::entity e) {
        auto item = old_to_new_map.find(e);

        if (item == old_to_new_map.end()) {
            auto new_e        = new_registry.create();
            old_to_new_map[e] = new_e;
            return new_e;
        }

        return item->second;
    }
};

template <class T>
void copy_marker_component(entt::registry const& from,
                           EntityMapper&         mapper,
                           entt::registry&       to) {
    for (auto [e] : from.view<T>().each()) {
        to.emplace_or_replace<T>(mapper(e));
    }
}

template <class T>
void copy_plain_component(entt::registry const& from,
                          EntityMapper&         mapper,
                          entt::registry&       to) {
    for (auto [e, c] : from.view<T>().each()) {
        to.emplace_or_replace<T>(mapper(e), c);
    }
}

template <class T>
void copy_nested_component(entt::registry const& from,
                           EntityMapper&         mapper,
                           entt::registry&       to) {
    for (auto [e, c] : from.view<T>().each()) {
        T component_copy = c;
        component_copy.remap_entities(mapper);
        to.emplace_or_replace<T>(mapper(e), component_copy);
    }
}

DatabaseCloneResult clone_database_with_entity_map(Database const& from,
                                                   QString  new_database_name,
                                                   QObject* p) {
    auto  ret           = std::make_unique<Database>(new_database_name, p);
    auto& from_registry = from.as_registry();
    auto& to_registry   = ret->as_registry();

    EntityMapper mapper {
        .new_registry   = to_registry,
        .old_to_new_map = { },
    };

    copy_marker_component<InvisibleComponent>(
        from_registry, mapper, to_registry);

    copy_marker_component<DisabledComponent>(
        from_registry, mapper, to_registry);

    copy_marker_component<VirtualTagComponent>(
        from_registry, mapper, to_registry);

    copy_plain_component<IdentityComponent>(from_registry, mapper, to_registry);


    copy_marker_component<ElementComponent>(from_registry, mapper, to_registry);


    copy_nested_component<ChildOfComponent>(from_registry, mapper, to_registry);

    copy_nested_component<ChildrenComponent>(
        from_registry, mapper, to_registry);


    copy_plain_component<TransformComponent>(
        from_registry, mapper, to_registry);

    copy_plain_component<GlobalTransformComponent>(
        from_registry, mapper, to_registry);


    ret->ray_source_resource.set(from.ray_source_resource.require().clone());
    ret->simulation_parameters_resource.set(
        from.simulation_parameters_resource.require());


    copy_plain_component<MaterialComponent>(from_registry, mapper, to_registry);

    copy_nested_component<MaterialGroupComponent>(
        from_registry, mapper, to_registry);

    copy_nested_component<MaterialGroupMemberComponent>(
        from_registry, mapper, to_registry);

    // GeometryComponent
    {
        for (auto [e, c] : from_registry.view<GeometryComponent>().each()) {
            auto local = c.clone();
            to_registry.emplace<GeometryComponent>(mapper(e), local);
        }
    }

    copy_nested_component<GeometryGroupComponent>(
        from_registry, mapper, to_registry);

    copy_nested_component<GeometryGroupMemberComponent>(
        from_registry, mapper, to_registry);

    copy_marker_component<TagComponent>(from_registry, mapper, to_registry);

    // ATagMemberComponent
    {
        // we need to know what all tags are out there

        auto all_tags = QSet<entt::entity>();

        for (auto [e, c] :
             from_registry.view<TagMembershipComponent>().each()) {
            for (auto const& t : c.tags) {
                all_tags.insert(t);
            }
        }

        for (auto tag_ent : std::as_const(all_tags)) {
            auto from_storage = from_registry.storage<ATagMemberComponent>(
                entt::to_integral(tag_ent));

            auto& to_storage = to_registry.storage<ATagMemberComponent>(
                entt::to_integral(mapper(tag_ent)));

            for (auto element : *from_storage) {
                to_storage.emplace(mapper(element));
            }
        }
    }

    copy_nested_component<TagMembershipComponent>(
        from_registry, mapper, to_registry);

    copy_plain_component<ImportErrorComponent>(
        from_registry, mapper, to_registry);

    copy_marker_component<SelectedComponent>(
        from_registry, mapper, to_registry);

    copy_plain_component<ColorComponent>(from_registry, mapper, to_registry);

    copy_plain_component<HasFluxMapComponent>(
        from_registry, mapper, to_registry);

    return DatabaseCloneResult {
        .database       = std::move(ret),
        .old_to_new_map = std::move(mapper.old_to_new_map),
    };
}

} // namespace db