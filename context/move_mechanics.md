# Move Mechanics

This document explains how polymer moves work: bond addition, kink, and removal.

## Move Pipeline

A move follows this pipeline in `move!(p::Polymer)`:

1. **Sample bond**: `p.moveIdx = rand(1:p.length - 1)` selects a chain position
2. **Identify bond**: `identifyBond!(p)` populates keys and coordinates
3. **Propose move**: `proposeMove!(p)` samples move type and new positions
4. **Check constraints**: overlap check, max length constraint
5. **Accept/reject**: Metropolis criterion via `acceptMove(p)`
6. **Execute**: `makeMove!(p)` calls appropriate function

## `identifyBond!(p)` Details

This function sets up pointers for the selected bond (see `move.jl:125-150`):

```julia
p.ki = p.keys_vec[p.moveIdx]     # key of bead i
p.kj = p.pos[p.ki][5]            # key of bead j (next bead)
p.kprev = p.pos[p.ki][4]         # key before bead i
p.knext = p.kj == 0 ? 0 : p.pos[p.kj][5]  # key after bead j

p.bi = view(p.pos[p.ki], 1:3)    # coordinates of bead i
p.bj = ...                       # coordinates of bead j
p.bprev = ...                    # coordinates of previous bead
p.bnext = ...                    # coordinates of next bead
```

## Move Types

| `moveType` | Name | Description |
|------------|------|-------------|
| 1 | Bond Addition | Insert two new beads between existing bond |
| 2 | Kink | Move one bead perpendicular to bond direction |
| 3 | Bond Removal | Remove two beads, joining neighbors directly |

### Bond Addition (`addBond!`)

Adds two new beads at `newPosi` and `newPosj`. The chain grows by 2 beads:
```
before: ... - ki - kj - ...
after:  ... - ki - new1 - new2 - kj - ...
```

Key operations:
- Allocate new keys via `KEYCOUNT`
- Create new entries in `pos` with proper prev/next pointers
- Re-route existing pointers: `p.pos[p.ki][5] = next_key1` and `p.pos[p.kj][4] = next_key2`
- Push new keys to `keys_vec`

### Bond Removal (`removeBond!`)

Removes beads at `ki` and `kj`, joining their neighbors:
```
before: ... - kprev - ki - kj - knext - ...
after:  ... - kprev - knext - ...
```

Key operations:
- Set lattice sites to 0
- Re-route pointers: `p.pos[p.kprev][5] = p.knext` and `p.pos[p.knext][4] = p.kprev`
- Handle head update if removing from start
- Delete keys from `pos` and filter from `keys_vec`

### Kink (`kinkBond!`)

Moves one bead without changing chain length. Either bead i moves to `newPosj` or bead j moves to `newPosi`.

## Rejection Conditions

In `move!(p)`, a move is rejected if:
- `p.kj == 0` (selected the tail—no bond to work with)
- `overlap(p)` returns true (new positions conflict with existing beads)
- `p.moveType == 1 && p.length == p.maxlength` (can't add at max length)
- Metropolis criterion fails
