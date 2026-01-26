# Hamiltonian Calculations

## Full Hamiltonian: `ham!(p::Polymer)`

Located in `hamiltonian.jl:7-50`. Computes total energy by traversing the chain:

```julia
key = p.head
for i in 1:p.length
    pos = p.pos[key][1:3]
    
    # Chemical potential: -μ per bead
    totalChemEnergy -= m.chemPot
    
    # Interaction energy via nnp() neighbor counting
    for (sp, idx) in m.intSpecies
        nnp_ = nnp(m, sp, pos)  # count neighbors
        if self-interaction:
            totalIntEnergy += intMat[idx,idx] * nnp_ / 2.0  # avoid double-counting
        else:
            totalIntEnergy += intMat[m.intIdx, idx] * nnp_
        end
    end
    
    key = p.pos[key][5]  # next bead
end
```

Currently, this assumes all monomers are identical (`m = p.monomers[1]`). For heterogeneous polymers, the monomer type must be retrieved per-bead.

## Delta Hamiltonian: `deltaHam!(p::Polymer)`

Located in `hamiltonian.jl:53-185`. Computes energy change for a proposed move.

**Important**: This function reads `m = getPolyMonomer(p, p.moveIdx)` to get the monomer at the move location. The current implementation uses `moveIdx` (chain position) which indexes into `monomerIndex`.

### Bond Addition (moveType == 1)
- Chemical: `-2 * chemPot` (adding 2 beads)
- Interaction: count new neighbors at `newPosi` and `newPosj`, plus the new bond between them

### Kink (moveType == 2)
- Chemical: no change (same number of beads)
- Interaction: `(neighbors at new position) - (neighbors at old position)`

### Bond Removal (moveType == 3)
- Chemical: `+2 * chemPot` (removing 2 beads)
- Interaction: subtract neighbors at `bi` and `bj`, minus 1 (the bond between them)

## Cross-Species Interactions

When a polymer interacts with bulk species, the pattern involves:
1. Iterating over `m.intSpecies` tuples of `(species, canonical_index)`
2. Using `m.intMat[m.intIdx, idx]` for interaction strength
3. Calling `setDeltaHP(sp, val)` to propagate energy changes to partner species

## The `nnp()` Function

Counts nearest neighbors belonging to a given species at a position. Defined in `utils.jl`. This is how interaction energies are computed—by counting adjacent occupied sites.
