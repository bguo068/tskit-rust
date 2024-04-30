impl_tskteardown!(
    super::bindings::tsk_table_collection_t,
    super::bindings::tsk_table_collection_free
);
impl_tskteardown!(super::bindings::tsk_tree_t, super::bindings::tsk_tree_free);

impl_tskteardown!(
    super::bindings::tsk_edge_table_t,
    super::bindings::tsk_edge_table_free
);
impl_tskteardown!(
    super::bindings::tsk_node_table_t,
    super::bindings::tsk_node_table_free
);
impl_tskteardown!(
    super::bindings::tsk_site_table_t,
    super::bindings::tsk_site_table_free
);
impl_tskteardown!(
    super::bindings::tsk_mutation_table_t,
    super::bindings::tsk_mutation_table_free
);
impl_tskteardown!(
    super::bindings::tsk_individual_table_t,
    super::bindings::tsk_individual_table_free
);
impl_tskteardown!(
    super::bindings::tsk_population_table_t,
    super::bindings::tsk_population_table_free
);
impl_tskteardown!(
    super::bindings::tsk_provenance_table_t,
    super::bindings::tsk_provenance_table_free
);
impl_tskteardown!(
    super::bindings::tsk_migration_table_t,
    super::bindings::tsk_migration_table_free
);
impl_tskteardown!(
    super::bindings::tsk_treeseq_t,
    super::bindings::tsk_treeseq_free
);
