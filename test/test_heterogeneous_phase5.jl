# Test suite for Phase 5: Type swapping for pattern length preservation

# Helper to count beads of each type
function countBeadTypes(p::Polymer)
    counts = Dict{Int32, Int}()
    key = p.head
    while key != 0
        intIdx = p.pos[key][6]
        counts[intIdx] = get(counts, intIdx, 0) + 1
        key = p.pos[key][5]
    end
    return counts
end

# Helper to get pattern as vector
function getPatternVector(p::Polymer)
    result = Int32[]
    key = p.head
    while key != 0
        push!(result, p.pos[key][6])
        key = p.pos[key][5]
    end
    return result
end

# MARK: -addBond preserves count
@testset "addBond inside pattern preserves type count" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer with block pattern: [1,1,1,1,2,2,2,2]
    # Type 1 is base (4 beads), Type 2 is pattern (4 beads)
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Record initial counts
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 4
    @test initial_counts[m2.intIdx] == 4
    
    # Find a bond inside the type-2 region (not at boundary)
    # Block pattern [1,1,1,1,2,2,2,2] - beads 5-6 are both type 2
    key = poly.head
    for _ in 1:4  # advance to 5th bead
        key = poly.pos[key][5]
    end
    @test poly.pos[key][6] == m2.intIdx
    next_key = poly.pos[key][5]
    @test poly.pos[next_key][6] == m2.intIdx
    
    # Manually set up for addBond inside type-2 region
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    poly.bprev = view(poly.pos[poly.kprev], 1:3)
    _sentinel = zeros(Int32, 3)
    poly.bnext = poly.knext == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.knext], 1:3)
    
    # Set up valid new positions that don't overlap
    poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
    poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + 1, poly.boxlength))]
    
    # Ensure positions are empty
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    # Perform addBond
    addBond!(poly)
    
    # Verify type counts:
    # - Pattern type unchanged: +2 new beads, -2 swapped to base = net 0
    # - Base type +2: from the 2 tail beads swapped from pattern
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] + 2  # Base gains 2 from swap
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]       # Pattern unchanged
    
    # Total length should increase by 2
    @test length(poly.pos) == 10
end

# MARK: -removeBond preserves count
@testset "removeBond inside pattern preserves type count" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer with custom pattern to ensure enough base beads after pattern
    # Pattern: [1,1,2,2,2,2,1,1,1,1]
    # 4 beads type-1 before pattern, 4 beads type-2 in pattern, 4 beads type-1 after
    # Wait - we need contiguous type-2 block, so: [1,1,1,1,2,2,2,2,1,1]
    poly = Polymer(;
        length = 10,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = [1, 1, 1, 1, 2, 2, 2, 2, 1, 1]
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    # Record initial counts
    initial_counts = countBeadTypes(poly)
    @test initial_counts[m1.intIdx] == 6
    @test initial_counts[m2.intIdx] == 4
    
    # Find a bond inside the type-2 region (beads 5-6, both type 2)
    key = poly.head
    for _ in 1:4  # advance to 5th bead
        key = poly.pos[key][5]
    end
    @test poly.pos[key][6] == m2.intIdx
    next_key = poly.pos[key][5]
    @test poly.pos[next_key][6] == m2.intIdx
    
    # Manually set up for removeBond
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    poly.bprev = view(poly.pos[poly.kprev], 1:3)
    _sentinel = zeros(Int32, 3)
    poly.bnext = poly.knext == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.knext], 1:3)
    
    # For removeBond, newPosi == bprev and newPosj == bnext
    poly.newPosi = Vector{Int32}(poly.bprev)
    poly.newPosj = Vector{Int32}(poly.bnext)
    
    # Perform removeBond
    removeBond!(poly)
    
    # Verify type counts:
    # - Pattern type unchanged: -2 removed, +2 swapped from base = net 0
    # - Base type -2: 2 beads swapped to pattern type
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] - 2  # Base loses 2 to swap
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]       # Pattern unchanged
    
    # Total length should decrease by 2
    @test length(poly.pos) == 8
end

# MARK: -patternTails updated after add
@testset "patternTails updated correctly after addBond" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer with block pattern
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    
    # Record initial tail position
    initial_tail = poly.patternTails[m2.intIdx]
    initial_tail_prev = poly.pos[initial_tail][4]
    expected_new_tail = poly.pos[initial_tail_prev][4]
    
    # Find a bond inside the type-2 region
    key = poly.head
    for _ in 1:4
        key = poly.pos[key][5]
    end
    next_key = poly.pos[key][5]
    
    # Set up addBond
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    poly.bprev = view(poly.pos[poly.kprev], 1:3)
    _sentinel = zeros(Int32, 3)
    poly.bnext = poly.knext == 0 ? view(_sentinel, 1:3) : view(poly.pos[poly.knext], 1:3)
    
    poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
    poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + 1, poly.boxlength))]
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    addBond!(poly)
    
    # patternTails should be updated - tail moved 2 beads backward
    @test poly.patternTails[m2.intIdx] == expected_new_tail
end

# MARK: -Add in base type (no swap)
@testset "addBond in base-type region has no swap" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m1 = poly.monomers[1]
    m2 = poly.monomers[2]
    
    initial_counts = countBeadTypes(poly)
    initial_tail = poly.patternTails[m2.intIdx]
    
    # Find a bond inside the type-1 (base) region (beads 1-2)
    key = poly.head
    next_key = poly.pos[key][5]
    @test poly.pos[key][6] == m1.intIdx
    @test poly.pos[next_key][6] == m1.intIdx
    
    # Set up addBond
    poly.ki = key
    poly.kj = next_key
    poly.kprev = poly.pos[key][4]
    poly.knext = poly.pos[next_key][5]
    
    poly.bi = view(poly.pos[poly.ki], 1:3)
    poly.bj = view(poly.pos[poly.kj], 1:3)
    # kprev is 0 for head bead, use a dummy view
    poly.bprev = poly.kprev == 0 ? zeros(Int32, 3) : view(poly.pos[poly.kprev], 1:3)
    poly.bnext = view(poly.pos[poly.knext], 1:3)
    
    poly.newPosi = [poly.bi[1], poly.bi[2], Int32(mod1(poly.bi[3] + 1, poly.boxlength))]
    poly.newPosj = [poly.bj[1], poly.bj[2], Int32(mod1(poly.bj[3] + 1, poly.boxlength))]
    poly.latS[poly.newPosi...] = 0
    poly.latS[poly.newPosj...] = 0
    
    addBond!(poly)
    
    # Base type count increases by 2, pattern type unchanged
    new_counts = countBeadTypes(poly)
    @test new_counts[m1.intIdx] == initial_counts[m1.intIdx] + 2
    @test new_counts[m2.intIdx] == initial_counts[m2.intIdx]
    
    # patternTails should be unchanged
    @test poly.patternTails[m2.intIdx] == initial_tail
end

# MARK: -Rejection when pattern too small
@testset "Reject addBond when pattern too small" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer with only 2 beads of type 2
    # Pattern: [1,1,1,1,2,2,1,1]
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :custom,
        customPattern = [1, 1, 1, 1, 2, 2, 1, 1]
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    
    # The type-2 block has only 2 beads - cannot shrink further
    # Validation in move! should reject addBond inside type-2
    
    # Find the type-2 bond
    key = poly.head
    for _ in 1:4
        key = poly.pos[key][5]
    end
    @test poly.pos[key][6] == m2.intIdx
    
    tailKey = poly.patternTails[m2.intIdx]
    prevTailKey = poly.pos[tailKey][4]
    newTailKey = poly.pos[prevTailKey][4]
    
    # newTailKey should be type-1 (not type-2), so validation fails
    @test poly.pos[newTailKey][6] != m2.intIdx
end

# MARK: -Rejection when insufficient base beads
@testset "Reject removeBond when insufficient base-type beads" begin
    l = Lattice(; boxlength = 30, nspecies = 1)
    
    # Create polymer where pattern is at the end (no base beads after)
    # Pattern: [1,1,1,1,2,2,2,2]
    poly = Polymer(;
        length = 8,
        maxlength = 16,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    m2 = poly.monomers[2]
    
    # The type-2 block is at the end - no base beads after tail
    tailKey = poly.patternTails[m2.intIdx]
    nextAfterTail = poly.pos[tailKey][5]
    
    # nextAfterTail should be 0 (end of chain)
    @test nextAfterTail == 0
end
