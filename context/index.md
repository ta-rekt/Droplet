# Droplet Context Index

This directory contains context files to help agents and humans understand the codebase and implement features. Keep this index up to date as context evolves.

## Files

- [polymer_data_structures.md](polymer_data_structures.md) - How polymers store bead positions, keys, and chain traversal logic
- [move_mechanics.md](move_mechanics.md) - Bond addition, removal, and kink moves
- [hamiltonian_calculations.md](hamiltonian_calculations.md) - Energy calculation patterns for polymers and monomers

## Code Pointers

| Concept | File | Key Functions/Structs |
|---------|------|----------------------|
| Polymer type | `src/types.jl:73-135` | `Polymer` struct |
| PolyMonomer type | `src/types.jl:48-70` | `PolyMonomer` struct |
| Chain traversal | `src/hamiltonian.jl:7-48` | `ham!` function |
| Bond operations | `src/move.jl:15-120` | `addBond!`, `kinkBond!`, `removeBond!` |
| Monomer initialization | `src/init.jl:135-230` | `initMonomers!` function |
| Get monomer by index | `src/utils.jl:102-108` | `getPolyMonomer` function |
