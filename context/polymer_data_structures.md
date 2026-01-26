# Polymer Data Structures

## The `pos` Dictionary

The polymer's bead positions are stored in an `OrderedDict{Int32, Vector{Int32}}` called `pos`. Each entry maps a **key** (bead identifier) to a **5-element vector**:

```julia
pos[key] = [x, y, z, prev_key, next_key]
```

- Elements 1-3: xyz coordinates of the bead
- Element 4: key of the previous bead in the chain (0 if this is the head)
- Element 5: key of the next bead in the chain (0 if this is the tail)

Keys are not contiguous—they are assigned incrementally by `KEYCOUNT` and never reused after removal. This means you **cannot iterate over the dictionary directly** to traverse the chain in order. Instead, you must follow the linked-list structure starting from `p.head`.

## Chain Traversal Pattern

To visit beads in sequence order, use the following pattern (from `ham!` in `hamiltonian.jl`):

```julia
key = p.head                      # start at head
for i in 1:p.length
    pos = p.pos[key][1:3]         # get coordinates
    # ... do work with this bead ...
    key = p.pos[key][5]           # follow to next bead
end
```

The traversal ends when `key == 0` (reached the tail's next pointer).

## The `monomerIndex` Vector

Currently, monomer types are stored in a separate `monomerIndex::Vector{Int32}` of length `maxlength`. The index into this vector corresponds to the **position in `keys_vec`**, not the key itself.

To get the monomer type for a bead at chain position `idx`:
```julia
p.monomers[p.monomerIndex[idx]]  # see getPolyMonomer() in utils.jl
```

## The `keys_vec` Vector

`keys_vec::Vector{Int32}` stores all active keys. It's maintained during add/remove operations:
- `addBond!`: pushes new keys
- `removeBond!`: filters out removed keys

The vector is used to map `moveIdx` (chain position) to the actual key:
```julia
p.ki = p.keys_vec[p.moveIdx]
```

## Head and Tail

- `p.head`: key of the first bead
- Tail: the bead where `p.pos[key][5] == 0`

Both ends can change during bond removal at the ends of the chain.
