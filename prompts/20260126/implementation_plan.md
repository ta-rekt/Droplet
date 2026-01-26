# Implementation Plan: Heterogeneous Polymers

## Overview

This plan implements heterogeneous polymers with multiple monomer types, where patterns of specific types maintain fixed length while the base type region can grow/shrink.

---

## Phase 1: Store Monomer Type in `pos` (6-element vectors)

### 1.1 Modify `pos` vector structure

**Goal**: Change `pos` entries from 5-element to 6-element vectors, where element 6 is the monomer type.

**Files to modify**:
- `src/types.jl`: Update comments/documentation for `pos` field
- `src/init.jl`: All functions that create `pos` entries

**Changes**:

1. **`initUniform!`** (init.jl:17-31): Change `zeros(5)` → `zeros(6)`, set `pos[i][6] = 1`
2. **`initLine!`** (init.jl:33-48): Same pattern
3. **`initBox!`** (init.jl:50-90): Same pattern
4. **`initShape!(::Polymer)`** (init.jl:235-245): Already calls `initLine!`, no direct changes needed

5. **`addBond!`** (move.jl:15-57):
   - Change `Vector{Int32}(undef, 5)` → `Vector{Int32}(undef, 6)`
   - Set `p.pos[next_key1][6]` and `p.pos[next_key2][6]` to match the type of the current bond's beads (inherit from parent)

6. **`kinkBond!`** (move.jl:59-90):
   - Change `Vector{Int32}(undef, 5)` → `Vector{Int32}(undef, 6)`
   - Preserve monomer type: copy element 6 from old entry

7. **`removeBond!`** (move.jl:92-120): No changes needed—just deletes entries

### 1.2 Remove `monomerIndex` vector

**Goal**: Eliminate `monomerIndex` since type is now stored in `pos`.

**Files to modify**:
- `src/types.jl`: Remove `monomerIndex` field from `Polymer` struct
- `src/init.jl`: Remove `monomerIndex` initialization from `initMonomers!`
- `src/utils.jl`: Update `getPolyMonomer` to read type from `pos`

**Changes to `getPolyMonomer`** (utils.jl:102-108):

```julia
function getPolyMonomer(p::Polymer, beadKey::Int32)
    if p.monomers === nothing
        return nothing
    end
    monomerType = p.pos[beadKey][6]
    return p.monomers[monomerType]
end
```

**Important**: This changes the semantics—`beadKey` is now the dictionary key, not the chain position. All callers must be updated:
- `deltaHam!` in hamiltonian.jl uses `p.moveIdx` → should use `p.ki` (the key)

### 1.3 Update Hamiltonian calculations

**`ham!`** (hamiltonian.jl:7-50): Currently uses `m = p.monomers[1]` assuming homogeneous. Change to:

```julia
key = p.head
for i in 1:p.length
    pos = p.pos[key][1:3]
    m = p.monomers[p.pos[key][6]]  # get monomer type from pos
    # ... rest unchanged but uses per-bead m ...
    key = p.pos[key][5]
end
```

**`deltaHam!`** (hamiltonian.jl:53-185): Update call to `getPolyMonomer`:

```julia
m = getPolyMonomer(p, p.ki)  # use key, not moveIdx
```

### 1.4 Tests for Phase 1

Create `test/test_heterogeneous_phase1.jl`:

1. **Test 6-element pos vectors**: Initialize polymer, verify all entries have length 6
2. **Test type storage**: Initialize with custom pattern, verify types stored correctly in `pos[key][6]`
3. **Test chain traversal with types**: Walk chain, collect types, verify sequence
4. **Test `getPolyMonomer` with keys**: Verify correct monomer returned for various keys
5. **Test Hamiltonian consistency**: Compare `ham!` before/after with heterogeneous types
6. **Test `deltaHam!` consistency**: Verify `initial + delta == final` for add/kink/remove

---

## Phase 2: Pattern Initialization with Chain Traversal

### 2.1 Update `initMonomers!` to use chain traversal

**Goal**: Assign monomer types by traversing the chain and writing to `pos[key][6]`.

**Current approach** (init.jl:163-229): Fills `monomerIndex` vector directly.

**New approach**: Must traverse chain via `pos[key][5]` to set types in order.

```julia
function initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)
    # Create monomers vector (unchanged)
    p.monomers = Vector{PolyMonomer}(undef, p.nMonomerTypes)
    # ... create PolyMonomer instances ...
    
    # Build pattern vector based on p.pattern
    patternVec = buildPatternVector(p)  # helper function
    
    # Traverse chain and assign types
    key = p.head
    for i in 1:p.length
        p.pos[key][6] = patternVec[i]
        key = p.pos[key][5]
    end
    
    return p
end
```

**Helper function** `buildPatternVector(p::Polymer)`:
- Returns `Vector{Int32}` of length `p.length` with monomer types
- Implements :uniform, :alternating, :block, :random, :custom patterns
- This is essentially the current pattern logic, just returning a vector

### 2.2 Validation at initialization

Add a check that each pattern type (2+) appears in exactly one contiguous block:

```julia
function validateSinglePatternBlocks(patternVec::Vector{Int32})
    for t in 2:maximum(patternVec)
        indices = findall(x -> x == t, patternVec)
        if !isempty(indices) && indices != collect(indices[1]:indices[end])
            error("Pattern type $t must appear in exactly one contiguous block")
        end
    end
end
```

### 2.3 Tests for Phase 2

1. **Test uniform pattern**: All beads should have type 1
2. **Test custom pattern**: `[1,1,1,2,2,2,1,1,1]` → verify positions
3. **Test promoter-enhancer pattern**: `[1,1,2,2,1,1,3,3,1,1]` → verify types
4. **Test validation rejects non-contiguous**: `[1,2,1,2]` for type 2 should error

---

## Phase 3: Pattern Boundary Tracking (`patternTails`)

### 3.1 Add `patternTails` to Polymer struct

**File**: `src/types.jl`

Add to Polymer:
```julia
patternTails::Union{Nothing, Dict{Int32, Int32}} = nothing  # patternType -> tail key
```

### 3.2 Initialize `patternTails` after pattern assignment

In `initMonomers!`, after assigning types, populate `patternTails`:

```julia
p.patternTails = Dict{Int32, Int32}()
key = p.head
for i in 1:p.length
    currentType = p.pos[key][6]
    nextKey = p.pos[key][5]
    
    if currentType >= 2  # pattern type
        # Check if this is the tail of this pattern
        if nextKey == 0 || p.pos[nextKey][6] != currentType
            p.patternTails[currentType] = key
        end
    end
    key = nextKey
end
```

### 3.3 Tests for Phase 3

1. **Test patternTails initialization**: For pattern `[1,1,2,2,2,1,1,3,3,1]`:
   - `patternTails[2]` should point to 3rd occurrence of type 2 (5th bead)
   - `patternTails[3]` should point to 2nd occurrence of type 3 (9th bead)
2. **Test patternTails with single pattern**: `[1,1,2,2,1,1]` → only `patternTails[2]` exists

---

## Phase 4: Reject Moves at Type Boundaries

### 4.1 Add boundary check in `move!`

**Goal**: If the selected bead and its next bead have different types, reject the move.

**File**: `src/move.jl`, in `move!(p::Polymer)` after `identifyBond!(p)`:

```julia
# Reject if bond spans type boundary
if p.kj != 0 && p.pos[p.ki][6] != p.pos[p.kj][6]
    return false
end
```

### 4.2 Tests for Phase 4

1. **Test boundary rejection**: With pattern `[1,1,2,2,1,1]`:
   - Select bond between beads 2-3 (types 1-2) → reject
   - Select bond between beads 4-5 (types 2-1) → reject
   - Select bond between beads 1-2 (types 1-1) → allow (if other conditions pass)
2. **Test statistics**: Run many moves, verify boundary bonds are never modified

---

## Phase 5: Type Swapping for Pattern Length Preservation

### 5.1 Modify `addBond!` to swap pattern tail types

**Logic**: When adding inside a pattern (type >= 2):
1. Find the last two beads of that pattern via `patternTails`
2. Change their type from X to 1
3. This compensates for the 2 new beads of type X being added

```julia
function addBond!(p::Polymer)
    currentType = p.pos[p.ki][6]
    
    # Existing add logic...
    # Set new beads to current type
    p.pos[next_key1][6] = currentType
    p.pos[next_key2][6] = currentType
    
    # If adding inside a pattern, swap tail types
    if currentType >= 2
        tailKey = p.patternTails[currentType]
        prevTailKey = p.pos[tailKey][4]
        
        # Swap last two beads of pattern to base type
        p.pos[tailKey][6] = 1
        p.pos[prevTailKey][6] = 1
        
        # Update patternTails to new tail
        # Walk backwards to find new tail
        newTail = p.pos[prevTailKey][4]
        while newTail != 0 && p.pos[newTail][6] != currentType
            newTail = p.pos[newTail][4]
        end
        if newTail != 0
            p.patternTails[currentType] = newTail
        else
            delete!(p.patternTails, currentType)  # pattern exhausted
        end
    end
    
    return p
end
```

**Wait**: This logic has an issue. After swapping, we need to find the *new* tail of the pattern. But since we're swapping at the tail, the new tail is actually two positions back. Let me reconsider...

Actually, the pattern shrinks by 2 at the tail (beads swapped to type 1) and grows by 2 at the location of bond insertion (new beads have type X). The net effect is that the pattern region "moves" but keeps the same count.

**Updating patternTails after add**:
```julia
# The new tail is the bead before prevTailKey with the same type
key = p.pos[prevTailKey][4]
while key != 0 && p.pos[key][6] == currentType
    key = p.pos[key][4]
end
# Now key is either 0 or has a different type; the previous key was the last of this type
# Actually we need to walk forward again... this is getting complex
```

**Simpler approach**: After the swap, re-scan to find the new tail:
```julia
if currentType >= 2
    # ... do the swap ...
    
    # Re-find tail by walking the pattern
    key = p.head
    newTail = nothing
    while key != 0
        if p.pos[key][6] == currentType
            newTail = key
        end
        key = p.pos[key][5]
    end
    if newTail !== nothing
        p.patternTails[currentType] = newTail
    else
        delete!(p.patternTails, currentType)
    end
end
```

This is O(n) but correct. Can be optimized later.

### 5.2 Modify `removeBond!` to swap base types into pattern

**Logic**: When removing inside a pattern (type >= 2):
1. Find the current tail and the bead after it
2. Change the bead after the tail (and the one after that) from type 1 to type X
3. This compensates for the 2 beads of type X being removed

```julia
function removeBond!(p::Polymer)
    currentType = p.pos[p.ki][6]
    
    # If removing inside a pattern, extend pattern at tail
    if currentType >= 2
        tailKey = p.patternTails[currentType]
        nextAfterTail = p.pos[tailKey][5]
        nextNextAfterTail = nextAfterTail == 0 ? 0 : p.pos[nextAfterTail][5]
        
        # Validate: there must be at least 2 base-type beads after the tail
        if nextAfterTail == 0 || p.pos[nextAfterTail][6] != 1 ||
           nextNextAfterTail == 0 || p.pos[nextNextAfterTail][6] != 1
            error("Cannot remove inside pattern: insufficient base-type beads at tail")
        end
        
        # Swap two base-type beads to pattern type
        p.pos[nextAfterTail][6] = currentType
        p.pos[nextNextAfterTail][6] = currentType
        
        # Update patternTails
        p.patternTails[currentType] = nextNextAfterTail
    end
    
    # Existing remove logic...
    # (delete entries for ki and kj)
    
    return p
end
```

### 5.3 Tests for Phase 5

1. **Test addBond inside pattern preserves count**: 
   - Start with `[1,1,2,2,2,2,1,1]` (4 type-2 beads)
   - Add bond inside type-2 region
   - Verify still exactly 4 type-2 beads
2. **Test removeBond inside pattern preserves count**:
   - Start with `[1,1,2,2,2,2,1,1]` (4 type-2 beads)
   - Remove bond inside type-2 region
   - Verify still exactly 4 type-2 beads
3. **Test patternTails updated correctly after add/remove**
4. **Test error on insufficient base-type beads for remove**

---

## Phase 6: Energy Accounting for Type Swaps

### 6.1 Extend `deltaHam!` to account for swap energy

**Goal**: When add/remove causes type swaps, include the energy change from those swaps.

For **addBond** inside pattern:
- Two beads at tail change from type X to type 1
- Energy change: loss of type-X interactions, gain of type-1 interactions

For **removeBond** inside pattern:
- Two beads after tail change from type 1 to type X
- Energy change: loss of type-1 interactions, gain of type-X interactions

**Implementation sketch** in `deltaHam!`:

```julia
if p.moveType == 1  # addBond
    currentType = p.pos[p.ki][6]
    if currentType >= 2
        # Find the two beads that will be swapped
        tailKey = p.patternTails[currentType]
        prevTailKey = p.pos[tailKey][4]
        tailPos = p.pos[tailKey][1:3]
        prevTailPos = p.pos[prevTailKey][1:3]
        
        mOld = p.monomers[currentType]  # monomer being swapped FROM
        mNew = p.monomers[1]            # monomer being swapped TO
        
        for (sp, idx) in mOld.intSpecies
            # Subtract old interactions
            p.deltaH -= mOld.intMat[mOld.intIdx, idx] * nnp(p, sp, tailPos)
            p.deltaH -= mOld.intMat[mOld.intIdx, idx] * nnp(p, sp, prevTailPos)
            
            # Add new interactions (with type-1 potentials)
            p.deltaH += mNew.intMat[mNew.intIdx, idx] * nnp(p, sp, tailPos)
            p.deltaH += mNew.intMat[mNew.intIdx, idx] * nnp(p, sp, prevTailPos)
        end
    end
end
```

Similar logic for `moveType == 3` (removeBond), but with opposite swap direction.

### 6.2 Tests for Phase 6

1. **Test deltaHam consistency with swaps**: 
   - Compute `ham!` before
   - Do add inside pattern (which causes swap)
   - Verify `initial + deltaH == ham!` after
2. **Test with different interaction potentials**: 
   - Set type-2 to have different `intPot` than type-1
   - Verify energy change is correct

---

## Summary Checklist

| Phase | Description | Files Modified | Tests |
|-------|-------------|----------------|-------|
| 1 | 6-element pos vectors, remove monomerIndex | types.jl, init.jl, move.jl, utils.jl, hamiltonian.jl | test_heterogeneous_phase1.jl |
| 2 | Pattern init with chain traversal | init.jl | test_heterogeneous_phase2.jl |
| 3 | patternTails tracking | types.jl, init.jl | test_heterogeneous_phase3.jl |
| 4 | Boundary rejection | move.jl | test_heterogeneous_phase4.jl |
| 5 | Type swapping in add/remove | move.jl | test_heterogeneous_phase5.jl |
| 6 | Energy accounting for swaps | hamiltonian.jl | test_heterogeneous_phase6.jl |

---

## Open Questions

1. **What happens when a pattern region exhausts base-type neighbors?** For example, if we keep removing from a pattern until there are no more type-1 beads to convert. Should this be an error, or should we prevent such moves?

2. **Multiple patterns of the same type**: The instructions say to check that there's only one pattern of any given type. But what about type 1 (base type)—is it okay for type-1 regions to be interspersed between patterns?

3. **Kink moves inside patterns**: The current plan only addresses add/remove for pattern length preservation. Do kink moves need special handling? (I assume not, since they don't change bead count.)

4. **`deltaHam!` for kink with heterogeneous types**: Currently deltaHam uses `getPolyMonomer(p, p.moveIdx)`. After the change to use keys, if a kink moves bead ki, should we use `p.pos[p.ki][6]` to get the type? (The type shouldn't change during a kink.)

5. **Chain position vs key confusion**: Currently `moveIdx` is a chain position (1 to length-1) used to index into `keys_vec`. After phase 1, `getPolyMonomer` will expect a key. Need to verify all call sites are updated correctly.

## Answers

1. Reject the move
2. Yes, type 1 is an exception. It can be interspersed between other types because we don't need to track its tail in the patternTail dictionary.
3. Yes, kinks don't need pattern-preservation logic.
4. Notice that you can find the new tail easily by moving backwards along the polymer. p.pos[key][4] contains the key of the previous bead, just like p.pos[key][5] points to the next bead. you can use this to update the keys in patternTails with the keys two beads prior for an addition, and two beads after for a removal. This is O(1) and doesn't add much complexity.
5. Yes, this is the intended behavior.

# General guidelines for code editing

- keep changes as frugal as you can, keeping close to the style and naming conventions of the existing code. make sure you depart as little as possible from what is necessary to implement the required changes, and avoid code bloat wherever possible. this code should still be readable and maintainable by humans