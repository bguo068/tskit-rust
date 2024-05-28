import msprime
import tskit

demography = msprime.Demography()
demography.add_population(name="A", initial_size=10_000)
demography.add_population(name="B", initial_size=5_000)
demography.add_population(name="C", initial_size=1_000)
# demography.set_migration_rate("A", "B", 0.01)
demography.add_population_split(time=1000, derived=["A", "B"], ancestral="C")
ts1: tskit.TreeSequence = msprime.sim_ancestry(
    samples={"A": 1, "B": 1},
    demography=demography,
    random_seed=12,
    recombination_rate=0.00001,
    sequence_length=300.0,
    # record_migrations=True
)
ts1 = msprime.sim_mutations(ts1, rate=0.0001, random_seed=12)

ts2 = ts1.keep_intervals([(10, 130)], record_provenance=False)
ts3 = ts1.keep_intervals([(10, 40), (100, 200)], record_provenance=False)


ts1.dump("testdata/1.trees")
ts2.dump("testdata/2.trees")
ts3.dump("testdata/3.trees")

