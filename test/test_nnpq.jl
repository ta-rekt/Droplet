# Test suite for nnpq function - Section 1 of three-way interaction implementation

using Test
using Droplet
using Random

# =============================================================================
# Helper functions
# =============================================================================

function setup_simple_lattice(boxlength=10)
    l = Lattice(; boxlength=boxlength, nspecies=3)
    return l
end

function setup_bulkmonomers(l, n1=10, n2=10)
    a = BulkMonomer(; name="A", length=n1)
    b = BulkMonomer(; name="B", length=n2)
    return a, b
end

# =============================================================================
# PART 1: Basic nnpq with BulkMonomers only
# =============================================================================

@testset "nnpq - BulkMonomer context (generic version)" begin
    l = setup_simple_lattice(10)
    a, b = setup_bulkmonomers(l, 5, 5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, -0.5)
    intPot!(b, c, -0.5)
    
    init!(l, [a, b, c])
    
    # Place particles manually for controlled test
    # Clear lattices first
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # Place A at position (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # Place B at neighbor (5,5,6)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 6
    b.latS[5, 5, 6] = 1
    
    # Place C at same position as A (5,5,5) - overlap
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    pos = Int32[5, 5, 5]
    
    # nnpq(a, b, c, pos) should count B neighbors when A and C overlap
    # A is at (5,5,5), C is at (5,5,5), B is at (5,5,6) which is neighbor
    result = nnpq(a, b, c, pos)
    @test result == 1  # B is a neighbor, and C overlaps at pos
    
    # Now test when C is NOT at the position
    c.latS[5, 5, 5] = 0
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 1, 1, 1
    c.latS[1, 1, 1] = 1
    
    result = nnpq(a, b, c, pos)
    @test result == 0  # C doesn't overlap, so coupling is zero
end

@testset "nnpq - multiple neighbors" begin
    l = setup_simple_lattice(10)
    a, b = setup_bulkmonomers(l, 5, 5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear lattices
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # Place A at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # Place C at same position (5,5,5) - overlap
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    # Place B at multiple neighbor positions
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 6  # +z neighbor
    b.latS[5, 5, 6] = 1
    b.pos[2][1], b.pos[2][2], b.pos[2][3] = 5, 6, 5  # +y neighbor
    b.latS[5, 6, 5] = 2
    b.pos[3][1], b.pos[3][2], b.pos[3][3] = 6, 5, 5  # +x neighbor
    b.latS[6, 5, 5] = 3
    
    pos = Int32[5, 5, 5]
    result = nnpq(a, b, c, pos)
    @test result == 3  # 3 B neighbors, C overlaps
end

# =============================================================================
# PART 2: nnpq with PolyMonomer as context
# =============================================================================

@testset "nnpq - PolyMonomer context checks intIdx" begin
    l = setup_simple_lattice(20)
    
    # Create polymer with two monomer types
    poly = Polymer(;
        length = 10,
        maxlength = 20,
        nMonomerTypes = 2,
        pattern = :uniform
    )
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    
    init!(l, [poly, a, b])
    
    # Get the two PolyMonomer types
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    
    # Manually set up a controlled configuration
    # Clear bulk lattices
    a.latS .= 0
    b.latS .= 0
    
    # Find a polymer bead of type 1
    pm1_x, pm1_y, pm1_z = Int32(0), Int32(0), Int32(0)
    for (key, posdata) in poly.pos
        if posdata[6] == pm1.intIdx
            pm1_x, pm1_y, pm1_z = posdata[1], posdata[2], posdata[3]
            break
        end
    end
    @test pm1_x > 0  # Should find at least one
    
    # Place A at same position as pm1 bead
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_x, pm1_y, pm1_z
    a.latS[pm1_x, pm1_y, pm1_z] = 1
    
    # Place B at a neighbor of pm1_pos
    neighbor_x = pm1_x + Int32(1)
    if neighbor_x > l.boxlength
        neighbor_x = Int32(1)
    end
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = neighbor_x, pm1_y, pm1_z
    b.latS[neighbor_x, pm1_y, pm1_z] = 1
    
    pm1_pos = Int32[pm1_x, pm1_y, pm1_z]
    
    # nnpq with pm1 as context should count because A overlaps with pm1
    result1 = nnpq(a, b, pm1, pm1_pos)
    @test result1 == 1
    
    # nnpq with pm2 as context should return 0 because A doesn't overlap pm2
    # (the lattice at pm1_pos has pm1.intIdx, not pm2.intIdx)
    result2 = nnpq(a, b, pm2, pm1_pos)
    @test result2 == 0
end

@testset "nnpq - heterogeneous polymer distinguishes types" begin
    l = setup_simple_lattice(30)
    
    # Create polymer with block pattern [1,1,...,1,2,2,...,2]
    poly = Polymer(;
        length = 10,
        maxlength = 20,
        nMonomerTypes = 2,
        pattern = :block
    )
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, poly.monomers[1], 0.0)
    intPot!(a, poly.monomers[2], 0.0)
    intPot!(b, poly.monomers[1], 0.0)
    intPot!(b, poly.monomers[2], 0.0)
    
    # Save monomer references before init! (which re-keys the dict)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    
    init!(l, [poly, a, b])
    
    # Clear bulk lattices
    a.latS .= 0
    b.latS .= 0
    
    # Find positions of both types by walking the chain
    pm1_x, pm1_y, pm1_z = Int32(0), Int32(0), Int32(0)
    pm2_x, pm2_y, pm2_z = Int32(0), Int32(0), Int32(0)
    
    key = poly.head
    while key != 0
        posdata = poly.pos[key]
        if posdata[6] == pm1.intIdx && pm1_x == 0
            pm1_x, pm1_y, pm1_z = posdata[1], posdata[2], posdata[3]
        elseif posdata[6] == pm2.intIdx && pm2_x == 0
            pm2_x, pm2_y, pm2_z = posdata[1], posdata[2], posdata[3]
        end
        key = posdata[5]
    end
    
    # Block pattern should have both types
    if pm1_x > 0 && pm2_x > 0
        # Place A at a pm1 position
        a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_x, pm1_y, pm1_z
        a.latS[pm1_x, pm1_y, pm1_z] = 1
        
        # Place B as neighbor
        neighbor_z = pm1_z + Int32(1)
        if neighbor_z > l.boxlength
            neighbor_z = Int32(1)
        end
        b.pos[1][1], b.pos[1][2], b.pos[1][3] = pm1_x, pm1_y, neighbor_z
        b.latS[pm1_x, pm1_y, neighbor_z] = 1
        
        pm1_pos = Int32[pm1_x, pm1_y, pm1_z]
        
        # Test: nnpq correctly identifies pm1 context
        @test nnpq(a, b, pm1, pm1_pos) == 1
        @test nnpq(a, b, pm2, pm1_pos) == 0  # pm2 not at this position
        
        # Now test with pm2 position
        a.latS[pm1_x, pm1_y, pm1_z] = 0
        a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm2_x, pm2_y, pm2_z
        a.latS[pm2_x, pm2_y, pm2_z] = 1
        
        neighbor_z2 = pm2_z + Int32(1)
        if neighbor_z2 > l.boxlength
            neighbor_z2 = Int32(1)
        end
        b.latS[pm1_x, pm1_y, neighbor_z] = 0
        b.pos[1][1], b.pos[1][2], b.pos[1][3] = pm2_x, pm2_y, neighbor_z2
        b.latS[pm2_x, pm2_y, neighbor_z2] = 1
        
        pm2_pos = Int32[pm2_x, pm2_y, pm2_z]
        
        @test nnpq(a, b, pm2, pm2_pos) == 1
        @test nnpq(a, b, pm1, pm2_pos) == 0
    else
        @warn "Block pattern didn't create both monomer types - skipping test"
        @test_skip true
    end
end

# =============================================================================
# PART 3: nnpq with PolyMonomer as partner
# =============================================================================

@testset "nnpq - PolyMonomer partner checks intIdx" begin
    l = setup_simple_lattice(20)
    
    poly = Polymer(;
        length = 10,
        maxlength = 20,
        nMonomerTypes = 2,
        pattern = :uniform
    )
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    
    a = BulkMonomer(; name="A", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(c, -1.0)
    
    init!(l, [poly, a, c])
    
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    
    # Clear bulk lattices
    a.latS .= 0
    c.latS .= 0
    
    # Find a pm1 bead position
    pm1_x, pm1_y, pm1_z = Int32(0), Int32(0), Int32(0)
    for (key, posdata) in poly.pos
        if posdata[6] == pm1.intIdx
            pm1_x, pm1_y, pm1_z = posdata[1], posdata[2], posdata[3]
            break
        end
    end
    @test pm1_x > 0
    
    # Place A at a position adjacent to pm1 bead
    test_x = pm1_x - Int32(1)
    if test_x < 1
        test_x = l.boxlength
    end
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = test_x, pm1_y, pm1_z
    a.latS[test_x, pm1_y, pm1_z] = 1
    
    # Place C at same position as A (overlap)
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = test_x, pm1_y, pm1_z
    c.latS[test_x, pm1_y, pm1_z] = 1
    
    test_pos = Int32[test_x, pm1_y, pm1_z]
    
    # nnpq(a, pm1, c, test_pos): pm1 is neighbor (at pm1_x,y,z), C overlaps
    result1 = nnpq(a, pm1, c, test_pos)
    @test result1 == 1
    
    # nnpq(a, pm2, c, test_pos): pm2 is NOT at pm1_pos, so should be 0
    result2 = nnpq(a, pm2, c, test_pos)
    @test result2 == 0
end

# =============================================================================
# PART 4: nnpq with both PolyMonomer partner and context
# =============================================================================

# Note: This test uses a block pattern where pm1 and pm2 beads are in separate blocks.
# The boundary between blocks will have adjacent beads of different types.
@testset "nnpq - both partner and context are PolyMonomer" begin
    l = setup_simple_lattice(30)
    
    # Block pattern: first half is type 1, second half is type 2
    poly = Polymer(;
        length = 10,
        maxlength = 20,
        nMonomerTypes = 2,
        pattern = :block
    )
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    
    a = BulkMonomer(; name="A", length=5)
    
    intPot!(a, -1.0)
    intPot!(a, poly.monomers[1], 0.0)
    intPot!(a, poly.monomers[2], 0.0)
    
    # Save monomer references before init! (which re-keys the dict)
    pm1 = poly.monomers[1]
    pm2 = poly.monomers[2]
    
    init!(l, [poly, a])
    
    # Clear bulk lattice
    a.latS .= 0
    
    # Find adjacent beads of different types (at block boundary)
    pm1_x, pm1_y, pm1_z = Int32(0), Int32(0), Int32(0)
    pm2_x, pm2_y, pm2_z = Int32(0), Int32(0), Int32(0)
    found = false
    
    key = poly.head
    while key != 0
        posdata = poly.pos[key]
        next_key = posdata[5]
        if next_key != 0
            next_posdata = poly.pos[next_key]
            if posdata[6] == pm1.intIdx && next_posdata[6] == pm2.intIdx
                pm1_x, pm1_y, pm1_z = posdata[1], posdata[2], posdata[3]
                pm2_x, pm2_y, pm2_z = next_posdata[1], next_posdata[2], next_posdata[3]
                found = true
                break
            end
        end
        key = posdata[5]
    end
    
    if found
        # Place A at the pm1 position (at the block boundary)
        a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_x, pm1_y, pm1_z
        a.latS[pm1_x, pm1_y, pm1_z] = 1
        
        pm1_pos = Int32[pm1_x, pm1_y, pm1_z]
        
        # nnpq(a, pm2, pm1, pm1_pos):
        # - partner pm2 is at adjacent position (next bead in chain)
        # - context pm1 is at pm1_pos (overlaps with A)
        result = nnpq(a, pm2, pm1, pm1_pos)
        @test result >= 1  # At least the adjacent pm2 bead
        
        # nnpq(a, pm1, pm2, pm1_pos):
        # - context pm2 is NOT at pm1_pos, so should be 0
        result2 = nnpq(a, pm1, pm2, pm1_pos)
        @test result2 == 0
    else
        @warn "Could not find adjacent beads of different types at block boundary"
        @test_skip true
    end
end

# =============================================================================
# PART 5: Boundary conditions
# =============================================================================

@testset "nnpq - periodic boundary conditions" begin
    l = setup_simple_lattice(10)
    a, b = setup_bulkmonomers(l, 5, 5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear lattices
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # Place A at boundary (1,1,1)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 1, 1, 1
    a.latS[1, 1, 1] = 1
    
    # Place C at same position
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 1, 1, 1
    c.latS[1, 1, 1] = 1
    
    # Place B at wrapped neighbor (boxlength, 1, 1) - wraps to neighbor of (1,1,1)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = l.boxlength, 1, 1
    b.latS[l.boxlength, 1, 1] = 1
    
    pos = Int32[1, 1, 1]
    result = nnpq(a, b, c, pos)
    @test result == 1  # B is neighbor across boundary, C overlaps
end

# =============================================================================
# PART 6: Edge cases
# =============================================================================

@testset "nnpq - no neighbors" begin
    l = setup_simple_lattice(10)
    a, b = setup_bulkmonomers(l, 5, 5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear lattices
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # Place A at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # Place C at same position
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    # Place B far away (no neighbors)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 1, 1, 1
    b.latS[1, 1, 1] = 1
    
    pos = Int32[5, 5, 5]
    result = nnpq(a, b, c, pos)
    @test result == 0  # No B neighbors
end

@testset "nnpq - no context overlap" begin
    l = setup_simple_lattice(10)
    a, b = setup_bulkmonomers(l, 5, 5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear lattices
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # Place A at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # Place B as neighbor
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 6
    b.latS[5, 5, 6] = 1
    
    # Place C elsewhere (no overlap with A)
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 1, 1, 1
    c.latS[1, 1, 1] = 1
    
    pos = Int32[5, 5, 5]
    result = nnpq(a, b, c, pos)
    @test result == 0  # C doesn't overlap, so coupling is zero
end

println("All nnpq tests completed!")
