SelectItemPopup {
    signal selectedElement(db_entity element)

    // We do this dance to avoid slowing down on large scenes
    // Only load content on open
    onOpened: source_model = AppData.layout.all_elements_model
    onClosed: source_model = null
    onSelectedEntity: (entity) => selectedElement(entity)
}
