# Integrated mechanistic interpretation

## Scope

This report integrates replica-reproducible protein distance changes with
ProLIF contacts, Z944 water bridges, MM/GBSA decomposition, ligand-to-pocket
distances, and ligand dihedrals. It intentionally does not interpret the pooled
PCA clusters as metastable states because neither global clustering nor
whole-replica classification generalized across independent simulations.

## Main result

Ligand binding does not produce one replica-independent global conformational
state in the current trajectories. Instead, it produces reproducible local and
inter-domain rearrangements embedded within substantial replica-dependent
global heterogeneity.

The strongest recurring structural feature is a contraction centered on Phe384
in Domain I toward the Domain-II pocket. For Z944 relative to apo, the leading
changes include:

- Phe384-Leu920: -1.29 Å mean difference
- Phe384-Phe917: -1.12 Å
- Phe384-Val916: -1.26 Å
- Phe384-Leu872: -0.97 Å
- Phe384-Phe875: -1.09 Å
- Phe384-Thr921: -1.02 Å

Every listed distance changes in the same direction in runs 1, 2, and 3.
mZ944 shows the same general Phe384/Domain-II contraction relative to apo,
although its wider determinant network includes stronger selectivity-filter
and Domain-III/IV rearrangements.

## Direct binding anchors

### Thr921

Thr921 has support from every available evidence layer:

- reproducible structural-distance changes in two system comparisons
- high ligand-contact occupancy, reaching approximately 88% in the bound-system
  comparison table
- favorable MM/GBSA contribution of approximately -2.29 kcal/mol for Z944 and
  -2.38 kcal/mol for mZ944
- mean Z944 water-bridge occupancy of approximately 28.5% across the three runs

Thr921 is therefore a strong conserved coupling point between direct ligand
recognition, hydration, and protein geometry.

### Lys1462

Lys1462 also has support from all four evidence layers:

- reproducible participation in the structural determinant network
- substantial ligand-contact occupancy
- favorable MM/GBSA contribution, especially for mZ944
- mean Z944 water-bridge occupancy of approximately 35%

The large between-replica variation in its mZ944 energetic contribution means
that its importance is strong but conformationally heterogeneous.

### Asn388

Asn388 combines reproducible structural involvement, approximately 41% Z944
contact occupancy, favorable Z944 MM/GBSA contribution, and a mean Z944
water-bridge occupancy above 5%. Its structural score is lower than that of
Phe384, but the convergence of independent evidence makes it a credible local
binding-network residue.

### Phe384

Phe384 is the dominant structural hub in the reproducible distance analysis.
It has moderate direct ligand contact and favorable energy but little
water-bridge occupancy. Its role is most consistent with hydrophobic/direct
packing coupled to the Domain-II pocket rather than hydration-mediated binding.

### Domain-II interaction core

Leu920, Thr921, Phe917, and Phe956 provide the clearest direct-interaction core:

- Leu920 and Thr921 have strong favorable energetic contributions.
- Phe917 combines reproducible structural contraction with contact and energy
  support.
- Phe956 is the most persistent ProLIF contact and strongest favorable
  per-residue MM/GBSA contributor, although it appears less frequently among
  the top protein-distance pairs.

This distinction is useful: Phe956 is a direct binding anchor, whereas Phe384
is a structural coupling hub.

## Structural-network and selectivity-filter determinants

Several leading reproducible residues have little direct ProLIF, MM/GBSA, or
water-bridge evidence. These should be interpreted as structural-network
readouts rather than ligand-contact residues.

Important examples include:

- Trp1781 in Domain IV
- Trp925 and Asn926 in Domain II
- Ser1461, Asp1463, and Trp1465 in Domain III
- Glu354 and Val357 in Domain I

For mZ944 relative to Z944, the leading recurring changes are:

- increased separation between residues 925/926 and 1461/1463
- decreased separation between residues 1776-1778 and 1464/1465

These changes indicate a reproducible redistribution of inter-domain
selectivity-filter geometry even though the pooled global landscapes and
cluster identities do not reproduce.

## Relationship to ligand position and conformation

The structural determinant pattern is consistent with the completed
ligand-position analyses:

- Z944 pocket-to-ligand distance: 3.01 ± 0.63 Å
- mZ944 pocket-to-ligand distance: 4.50 ± 1.29 Å
- Z944 filter-to-ligand distance: 11.14 ± 0.58 Å
- mZ944 filter-to-ligand distance: 10.23 ± 1.31 Å

mZ944 is therefore more displaced from the pocket center while, on average,
lying closer to the filter center. Its larger standard deviations agree with
the replica-dependent global geometry and the broader structural determinant
network.

The ligand-dihedral analysis provides an internal conformational correlate:

- chi3 shifts from approximately 121 degrees in Z944 to 100 degrees in mZ944
- chi2 is the most heterogeneous and replica-dependent torsion
- mZ944 has greater chi7 variability
- Z944 has greater chi4 and chi5 variability

These torsional differences offer a plausible local coordinate linking ligand
conformation to the distinct Domain-II/III/IV distance rearrangements, but a
frame-matched correlation analysis is required before claiming causality.

## Energetic interpretation

The three-replica overall MM/GBSA means are:

- Z944: -23.30 ± 1.36 kcal/mol
- mZ944: -21.01 ± 4.16 kcal/mol

Z944 is slightly more favorable and substantially more reproducible. This is
consistent with its shorter and less variable pocket distance. The energetic
difference should not be attributed to a single residue: the current
per-residue GBIS values are pairwise residue-ligand terms and are not strictly
additive.

## Defensible biological model

The combined evidence supports the following working model:

1. Both ligands couple Phe384 in Domain I to the Domain-II binding pocket.
2. A conserved direct-interaction core involving Phe917, Leu920, Thr921,
   Phe956, and nearby residues stabilizes ligand occupancy.
3. Thr921 and Lys1462 connect direct binding to water-mediated interaction
   networks.
4. mZ944 adopts a more displaced and heterogeneous pocket geometry and produces
   reproducible redistribution among Domain-II, Domain-III, and Domain-IV
   selectivity-filter distances.
5. These local and inter-domain effects recur across replicas, whereas a single
   global ligand-specific conformational basin does not.

## Evidence limits and next experiment

- Water-bridge data currently exist only for Z944.
- ProLIF, MM/GBSA, COM-distance, and ligand-dihedral descriptors are unavailable
  for apo where they require a ligand.
- The integrated rankings are based on replica-consistent direction and effect
  size, not independent-frame hypothesis tests.
- Frames are temporally correlated.
- State kinetics and MSM construction remain unjustified without reproducible
  state correspondence and lag-time validation.

The most informative next calculation is frame-matched analysis within the
ligand-bound systems: correlate the leading protein distances with ligand COM
distance, chi2/chi3/chi7, ProLIF contacts, and water-bridge occupancy separately
for each replica, then retain only relationships with consistent direction
across replicas.

## Frame-matched coupling results

The protein-distance checkpoints were matched to 600 COM and ligand-dihedral
observations per replica. Sixty distance/COM relationships retain the same
correlation sign across all three runs, but only one has a minimum absolute
Spearman correlation above 0.2:

- Z944 distance 351-356 versus ligand-to-filter distance:
  rho = 0.23, 0.34, and 0.22 in runs 1, 2, and 3.

The corresponding 351-356 distance is negatively related to Z944
ligand-to-pocket distance in all three runs (rho = -0.18, -0.29, and -0.24).
Thus, movement of Z944 toward the filter and away from the pocket is coupled to
a reproducible local Domain-I filter rearrangement. Most other COM
relationships are directionally reproducible but weak.

Circular-linear analysis reveals stronger torsional coupling:

- mZ944 chi2 couples to distances linking residues 1776-1778 with 1464/1465.
- mZ944 chi3 couples to the Phe384-1463/1465 structural network.
- Z944 chi2 couples to Phe384 distances involving residues 916, 918-920, and
  926.
- Z944 chi4 couples to Phe384-Thr921 and related Domain-II distances.
- mZ944 chi7 couples reproducibly to the 356-925 distance.

These results connect the previously observed ligand torsional differences to
specific protein/filter coordinates. They remain correlations and do not
establish whether ligand torsion drives protein motion or responds to it.

Z944 water-bridge coupling was evaluated for runs 2 and 3. The relationships
have consistent directions for many distances but are weak: the strongest
minimum absolute correlation is approximately 0.16 for the 956-1464 distance
versus bridge-event count. This supports a distributed hydration network rather
than a single dominant protein-distance gating coordinate. Run1 water
time-series data were not available in the same frame-resolved format and were
not imputed.

## Frame-matched ProLIF coupling

ProLIF interaction types were aggregated into residue-level any-contact
variables and matched to the distance checkpoints. Only relationships with
5-95% contact occupancy in every replica were retained. Of 27 testable
endpoint-contact relationships, seven have the same correlation sign in all
three runs, but none has a minimum absolute rho of 0.2.

The strongest recurring example is Z944 distance 387-955 versus Ile387 contact
presence (rho = -0.10, -0.19, and -0.12). Several Phe384-centered Z944
relationships are also consistently negative, but run1 correlations are close
to zero even where runs 2 and 3 are moderate.

Therefore, direct contact formation is not a robust one-dimensional switch
along any selected C-alpha distance. This does not contradict the strong
replica-averaged ProLIF occupancies or energetic contributions. Instead, it
indicates that contact formation depends on multivariate local geometry,
side-chain orientation, ligand pose, and hydration that are not captured by one
C-alpha distance.

The ProLIF fingerprints are sparse, sampled at approximately 3-6 ns, so this
analysis has lower temporal resolution than the COM and dihedral coupling
analysis.

## Structural visualization package

The top 25 reproducible pairs for each comparison have been mapped to
comparison-specific PDB structures. Signed, rank-weighted standardized effects
are stored in the B-factor field:

- negative/blue: net reproducible contraction
- positive/red: net reproducible expansion
- magnitude: normalized residue-level structural effect

The Z944-versus-apo mapping emphasizes the Phe384-centered contraction toward
Domain II. The mZ944-versus-Z944 network separates expanded Domain-II/III
connections involving residues 925/926 and 1461/1463 from contracted
Domain-III/IV connections involving residues 1464/1465 and 1776-1778.

Native VMD and PyMOL scripts display the mapped cartoon, determinant side
chains, ligand, and the ten leading pair distances. Domain-organized PNG/PDF
network figures provide a two-dimensional summary. The Z944-versus-apo VMD
script was executed in text mode and loaded all selections and graphics without
errors.

## Uncertainty-validated determinant subset

The top 25 determinants per comparison were evaluated using two complementary
uncertainty criteria:

1. 10-ns non-overlapping block-bootstrap intervals within every trajectory.
2. A 95% t interval across the three independent replica-level mean
   differences.

Dual support requires all three run-level block intervals to exclude zero in
the common direction and the across-replica interval to exclude zero.

Supported counts are:

- Z944 versus apo: 15 of 25
- mZ944 versus apo: 10 of 25
- mZ944 versus Z944: 12 of 25

The central Phe384 contraction is retained. Examples for Z944 versus apo are:

- Phe384-Leu920: -1.29 Å; replica interval [-1.91, -0.67]
- Phe384-Val916: -1.26 Å; [-2.24, -0.27]
- Phe384-Ile876: -1.13 Å; [-1.78, -0.48]
- Phe384-Phe917: -1.12 Å; [-1.95, -0.28]
- Phe384-Thr921: -1.02 Å; [-1.90, -0.15]

The mZ944-versus-Z944 filter-network result is also retained:

- 1776-1464: -1.36 Å; [-2.38, -0.33]
- 1777-1464: -1.34 Å; [-2.57, -0.11]
- 925-1463: +0.89 Å; [0.16, 1.61]
- 926-1463: +0.87 Å; [0.17, 1.56]
- 926-1459: +0.66 Å; [0.57, 0.74]

Only three independent replicas are available, so the across-replica
intervals have two degrees of freedom and remain distribution-sensitive. The
block bootstrap quantifies temporal uncertainty within trajectories; it does
not create additional independent biological replicas.
