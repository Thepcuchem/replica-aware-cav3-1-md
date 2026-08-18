# AI-driven MD framework: validated results report

## Study design

The framework analyzes apo, Z944, and mZ944 Cav3.1 systems using three
independent trajectories per system. Common structural analyses use the final
300 ns of each trajectory at 0.2-ns sampling:

- 1,501 frames per replica
- 13,509 frames total
- 2,016 common C-alpha distances
- binding-pocket and selectivity-filter residues shared across all systems

Ligand-specific ProLIF, COM-distance, dihedral, water-bridge, and MM/GBSA
descriptors are used as explanatory layers rather than imputed into apo.

## Validation sequence

1. A one-trajectory PCA baseline showed visually separated systems but was
   explicitly labeled exploratory.
2. Nine-replica RMSD classification gave mean held-out balanced accuracy
   0.381.
3. Nine-replica continuous distances gave frame-level balanced accuracy
   0.387 and
   trajectory accuracy 0.333.
4. The 8-angstrom contact map gave frame-level balanced accuracy
   0.258
   and trajectory accuracy
   0.222.
5. Independently discovered two-state clusters failed centroid recurrence:
   cross-replica matched states were farther apart than different states within
   a reference replica.

These controls reject a single replica-independent global ligand-state model
for the current trajectories. A pooled global MSM is therefore not supported.

## Reproducible molecular determinants

Feature-level comparison retained residue-pair distances with the same effect
direction in run1, run2, and run3. The strongest recurring result is contraction
of Phe384 in Domain I toward the Domain-II pocket in both ligand-bound systems.
Z944-versus-apo contractions include Phe384 with Val916, Phe917, Leu920,
Thr921, Ile876, Phe875, Asn926, and Leu872.

mZ944 differs from Z944 through a coordinated filter-network redistribution:

- expansion between Domain-II residues 925/926 and Domain-III residues
  1459/1463
- contraction between Domain-III residues 1464/1465 and Domain-IV residues
  1776-1778

## Uncertainty-supported subset

Dual uncertainty support requires every run's 10-ns block-bootstrap interval
and the across-replica 95% interval to exclude zero in the common direction.

- Z944 versus apo:
  15/25
- mZ944 versus apo:
  10/25
- mZ944 versus Z944:
  12/25

Only three independent replicas are available, so across-replica intervals have
two degrees of freedom.

## Mechanistic integration

The integrated evidence supports two residue classes:

- Direct/contact-coupled anchors: Thr921, Lys1462, Asn388, Phe384, Leu920,
  Phe917, Phe956, Ile387, Val1820, Gln1816, and Leu872.
- Structural-network determinants: Trp1781, Trp925, Asn926, Ser1461, Asp1463,
  Trp1465, Glu354, Val357, and related filter residues.

Thr921 and Lys1462 connect direct ligand recognition, favorable energy, and
Z944 water-mediated interactions. Phe956 is the strongest persistent direct
contact/energetic anchor, while Phe384 is the dominant protein-geometry hub.

## Frame-matched coupling

60 distance/COM
relationships retain their sign across all three replicas, but only
1 reaches minimum
absolute rho of 0.2. Z944 distance 351-356 tracks ligand movement toward the
filter and away from the pocket.

Circular-linear relationships connect:

- mZ944 chi2 to the 1776-1778/1464-1465 filter network
- mZ944 chi3 to the Phe384/1463-1465 network
- Z944 chi2 and chi4 to the Phe384/Domain-II pocket network

Frame-matched ProLIF analysis yields
7 sign-consistent relationships
from 27 testable pairs,
but 0 reaches minimum
absolute rho of 0.2. Direct contact switching is therefore multivariate rather
than controlled by one C-alpha distance.

## Defensible conclusion

Ligand binding does not impose one universal global conformational basin across
the current replicas. Instead, Z944 and mZ944 generate reproducible local and
inter-domain distance changes within broader replica-dependent landscapes.
Z944 is associated with a compact, reproducible Phe384-to-Domain-II pocket
geometry. mZ944 retains this coupling but adds a distinct redistribution of
Domain-II, Domain-III, and Domain-IV filter geometry, consistent with its more
displaced and heterogeneous ligand pose.

State kinetics, transition pathways, and a pooled MSM should remain deferred
until structurally recurrent states and lag-time convergence are demonstrated.
