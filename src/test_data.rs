/// mimic the c simulate function in tskit c api document
/// https://tskit.dev/tskit/docs/stable/c-api.html#basic-forwards-simulator
pub mod simulation {
    use core::panic;

    use crate::{
        metadata::{MetadataError, MetadataRoundtrip, PopulationMetadata},
        EdgeId, IndividualId, MutationId, NodeFlags, NodeId, PopulationId, Position,
        SimplificationOptions, SiteId, TableCollection, TableSortOptions, TreeSequence,
        TreeSequenceFlags, TskitError,
    };
    use rand::{rngs::StdRng, Rng, SeedableRng};

    struct MyMeta {
        inner: String,
    }
    impl From<String> for MyMeta {
        fn from(value: String) -> Self {
            MyMeta { inner: value }
        }
    }
    impl<'a> From<&'a str> for MyMeta {
        fn from(value: &'a str) -> Self {
            MyMeta {
                inner: value.to_owned(),
            }
        }
    }

    // helper structs, impls and functions

    impl MetadataRoundtrip for MyMeta {
        fn encode(&self) -> Result<Vec<u8>, MetadataError> {
            Ok(self.inner.as_bytes().to_owned())
        }
        fn decode(md: &[u8]) -> Result<Self, MetadataError>
        where
            Self: Sized,
        {
            Ok(MyMeta {
                inner: String::from_utf8(md.to_owned()).unwrap(),
            })
        }
    }

    impl PopulationMetadata for MyMeta {}

    fn add_pop(tables: &mut TableCollection, name: &str) -> PopulationId {
        tables
            .add_population_with_metadata(&MyMeta::from(name))
            .unwrap()
    }

    fn add_ind(
        tables: &mut TableCollection,
        parent1: (NodeId, NodeId),
        parent2: (NodeId, NodeId),
    ) -> IndividualId {
        let parent1_ind = tables.nodes().individual(parent1.0).unwrap();
        let parent2_ind = tables.nodes().individual(parent2.0).unwrap();
        let flags = 0u32;
        let loc_null = None;
        tables
            .add_individual(flags, loc_null, [parent1_ind, parent2_ind])
            .unwrap()
    }

    fn find_parent(
        rng: &mut StdRng,
        parents: &[(NodeId, NodeId)],
        child_pop: PopulationId,
    ) -> ((NodeId, NodeId), PopulationId) {
        assert_eq!(parents.len() % 2, 0);
        let (pop_anc, pop_1, pop_2) = (0, 1, 2);
        let child_pop: i32 = child_pop.into();

        let pop_size = parents.len();
        let mut parent_pop = child_pop;

        let is_migrant = (child_pop != pop_anc) && rng.gen_bool(0.01);
        if is_migrant {
            parent_pop = if child_pop == pop_1 { pop_2 } else { pop_1 };
        };
        let parent = match parent_pop {
            // pop_anc
            0 => parents[rng.gen_range(0..pop_size)],
            // pop_1
            1 => parents[rng.gen_range(0..(pop_size / 2))],
            // pop_2
            2 => parents[rng.gen_range((pop_size / 2)..pop_size)],
            _ => panic!("wrong population id encountered"),
        };
        (parent, parent_pop.into())
    }

    fn find_breakpoint(rng: &mut StdRng, seqlen: Position) -> Position {
        // avoid breaking as edges
        let seqlen = f64::from(seqlen).floor() as usize;
        let sel = rng.gen_range(1..seqlen) as f64;
        Position::from(sel)
    }

    fn add_node(
        tables: &mut TableCollection,
        is_sample: bool,
        time: usize,
        pop: PopulationId,
        ind: IndividualId,
    ) -> NodeId {
        tables
            .add_node(
                if is_sample {
                    NodeFlags::new_sample()
                } else {
                    NodeFlags::default()
                },
                time as f64,
                pop,
                ind,
            )
            .unwrap()
    }

    fn add_edge(
        tables: &mut TableCollection,
        start: impl Into<Position>,
        end: impl Into<Position>,
        parent_node: NodeId,
        child_node: NodeId,
    ) -> EdgeId {
        tables
            .add_edge(start, end, parent_node, child_node)
            .unwrap()
    }

    fn find_overlaps<P>(start: P, end: P, intervals: &Vec<(P, P)>, out: &mut Vec<(P, P)>)
    where
        P: Into<Position> + Copy + PartialOrd,
    {
        // assert intervals is sorted
        assert!(intervals.iter().all(|(a, b)| *a <= *b));
        assert!(intervals
            .iter()
            .zip(intervals.iter().skip(1))
            .all(|(p1, p2)| p1.1 <= p2.0));
        // clear out
        out.clear();

        for (m, n) in intervals {
            // no overlap
            if (*n <= start) || (end <= *m) {
                continue;
            }
            let new_start = if *m < start { start } else { *m };
            let new_end = if *n < end { *n } else { end };
            out.push((new_start, new_end));
        }
    }

    fn find_mutation_pos<P>(rng: &mut StdRng, s: P, e: P) -> usize
    where
        P: Into<Position>,
    {
        let s = f64::from(Into::<Position>::into(s)).ceil() as usize;
        let e = f64::from(Into::<Position>::into(e)).floor() as usize;
        rng.gen_range(s..e)
    }

    fn calc_derived_state(site_last_mutation_order: &[usize], mut_pos: usize) -> [u8; 1] {
        [b'a'
            + match site_last_mutation_order[mut_pos] + 1 {
                x if x > 45 => 45u8,
                x => x as u8,
            }]
    }

    /// simulate diplid individual with migration between two subpopulations
    ///
    /// Both full_trees and trucated_trees will be generated
    pub fn simulate_two_treesequences<P>(
        seqlen: P,
        pop_size: usize,
        start_time: usize,
        split_time: usize,
        intervals: &[(P, P)],
        seed: u64,
    ) -> Result<(TreeSequence, TreeSequence), TskitError>
    where
        P: Into<Position> + Copy + PartialOrd,
    {
        let rng = &mut StdRng::seed_from_u64(seed);
        let intervals: Vec<(Position, Position)> = intervals
            .iter()
            .map(|(a, b)| ((*a).into(), (*b).into()))
            .collect();
        assert!(split_time < start_time);
        assert_eq!(pop_size % 2, 0);
        // tables without truncation
        let mut tables = TableCollection::new(seqlen).unwrap();
        // expected tables after truncation
        // it is built following `tables` except for positions for edge table
        let mut tr_tbls = TableCollection::new(seqlen).unwrap();

        let mut buffer = Vec::new();

        // add pop
        let pop_anc = add_pop(&mut tables, "ancestor");
        let pop_1 = add_pop(&mut tables, "pop1");
        let pop_2 = add_pop(&mut tables, "pop2");

        add_pop(&mut tr_tbls, "ancestor");
        add_pop(&mut tr_tbls, "pop1");
        add_pop(&mut tr_tbls, "pop2");

        // state variables for site/mutation tables
        let num_sites = f64::from(seqlen.into()) as usize;
        let mut site_last_mutation_order = vec![0usize; num_sites];

        let mut site_last_mutation_tables = vec![MutationId::NULL; num_sites];
        let mut site_last_mutation_tr_tbls = vec![MutationId::NULL; num_sites];

        let mut site_id_map_tables = vec![SiteId::NULL; num_sites];
        let mut site_id_map_tr_tbls = vec![SiteId::NULL; num_sites];

        // base population
        let mut parents = Vec::<(NodeId, NodeId)>::with_capacity(pop_size);
        for _ in 0..pop_size {
            const FLAGS: u32 = 0;
            let loc_null = None;
            let parent_ind = tables.add_individual(FLAGS, loc_null, None).unwrap();
            tr_tbls.add_individual(FLAGS, loc_null, None).unwrap();

            let parent_id = (
                add_node(&mut tables, false, start_time, pop_anc, parent_ind),
                add_node(&mut tables, false, start_time, pop_anc, parent_ind),
            );
            parents.push(parent_id);
            //
            add_node(&mut tr_tbls, false, start_time, pop_anc, parent_ind);
            add_node(&mut tr_tbls, false, start_time, pop_anc, parent_ind);
        }

        // offspring population
        let mut children = Vec::<(NodeId, NodeId)>::with_capacity(pop_size);

        for t in (0..start_time).rev() {
            for i in 0..pop_size {
                // select breakpoints
                let breakpoint1 = find_breakpoint(rng, seqlen.into());
                let breakpoint2 = find_breakpoint(rng, seqlen.into());

                // find child pop
                let mut child_pop = pop_anc;
                if t > split_time {
                    child_pop = if i < pop_size / 2 { pop_1 } else { pop_2 }
                }

                // find parents
                let (parent1, _parent1_pop) = find_parent(rng, &parents, child_pop);
                let (parent2, _parent2_pop) = find_parent(rng, &parents, child_pop);

                // add individual
                let child_ind = add_ind(&mut tables, parent1, parent2);
                add_ind(&mut tr_tbls, parent1, parent2);

                // add nodes
                let is_sample = t == 0;
                let child_id = (
                    add_node(&mut tables, is_sample, t, child_pop, child_ind),
                    add_node(&mut tables, is_sample, t, child_pop, child_ind),
                );

                add_node(&mut tr_tbls, is_sample, t, child_pop, child_ind);
                add_node(&mut tr_tbls, is_sample, t, child_pop, child_ind);

                // add edges, sites & mutations to both tables and tr_tabls
                let mu = 0.01f64;
                for (s, e, p, c) in [
                    (0.0.into(), breakpoint1, parent1.0, child_id.0),
                    (breakpoint1, seqlen.into(), parent1.1, child_id.0),
                    (0.0.into(), breakpoint2, parent2.0, child_id.1),
                    (breakpoint2, seqlen.into(), parent2.1, child_id.1),
                ] {
                    add_edge(&mut tables, s, e, p, c);

                    let mut_pos = find_mutation_pos(rng, s, e);
                    let mut mut_prob = f64::from(e - s) * mu;
                    if mut_prob > 1.0 {
                        mut_prob = 1.0;
                    }
                    let to_add_mut: bool = rng.gen_bool(mut_prob);
                    let derived_state = &calc_derived_state(&site_last_mutation_order, mut_pos);
                    let t = t as f64;

                    if to_add_mut {
                        // add site
                        let site_not_exist = site_id_map_tables[mut_pos] == SiteId::NULL;
                        if site_not_exist {
                            site_id_map_tables[mut_pos] =
                                tables.add_site(mut_pos as f64, Some(&[b'a'])).unwrap();
                        }
                        // add mutation
                        let parent_mut = site_last_mutation_tables[mut_pos];
                        let site = site_id_map_tables[mut_pos];
                        let new_mutation = tables
                            .add_mutation(site, c, parent_mut, t, Some(derived_state))
                            .unwrap();

                        site_last_mutation_tables[mut_pos] = new_mutation;
                        site_last_mutation_order[mut_pos] += 1;
                    }

                    find_overlaps(s, e, &intervals, &mut buffer);
                    for (s_, e_) in buffer.iter() {
                        add_edge(&mut tr_tbls, *s_, *e_, p, c);
                        let mut_pos_f = mut_pos as f64;

                        if to_add_mut && (*s_ <= mut_pos_f) && (*e_ > mut_pos_f) {
                            // add site
                            let site_not_exist = site_id_map_tr_tbls[mut_pos] == SiteId::NULL;
                            if site_not_exist {
                                site_id_map_tr_tbls[mut_pos] =
                                    tr_tbls.add_site(mut_pos as f64, Some(&[b'a'])).unwrap();
                            }
                            // add mutation
                            let parent_mut = site_last_mutation_tr_tbls[mut_pos];
                            let site = site_id_map_tr_tbls[mut_pos];
                            let new_mutation = tr_tbls
                                .add_mutation(site, c, parent_mut, t, Some(derived_state))
                                .unwrap();
                            site_last_mutation_tr_tbls[mut_pos] = new_mutation;
                        }
                    }
                }

                // add edges for tr_tbls
                children.push(child_id);
            }
            // NOTE: avoid simplifcation so that both tables and tr_tables share the same ids

            // set children as parents and clear children
            std::mem::swap(&mut children, &mut parents);
            children.clear();
        }

        tables.full_sort(TableSortOptions::all()).unwrap();
        tr_tbls.full_sort(TableSortOptions::all()).unwrap();

        // simplify
        let mut samples = Vec::<NodeId>::with_capacity(pop_size * 2);
        parents
            .iter()
            .for_each(|p| samples.extend([p.0, p.1].iter()));

        tables
            .simplify(&samples, SimplificationOptions::default(), false)
            .unwrap();

        tr_tbls
            .simplify(&samples, SimplificationOptions::default(), false)
            .unwrap();

        // build indices
        tables.build_index().unwrap();
        tr_tbls.build_index().unwrap();

        // to tree sequences
        let full_trees = TreeSequence::new(tables, TreeSequenceFlags::default()).unwrap();
        let truncated_trees = TreeSequence::new(tr_tbls, TreeSequenceFlags::default()).unwrap();

        Ok((full_trees, truncated_trees))
    }
}
