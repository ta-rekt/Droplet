# Test file for Section 4: ham! with three-way coupling for BulkMonomer
# Tests the full Hamiltonian calculation including coupling contributions

using Test
using Droplet

println("Testing ham! with three-way coupling (Section 4)")
println("=" ^ 50)

# Helper to set up a simple lattice
function setup_lattice(boxlength=20)
    return Lattice(; boxlength=boxlength, nspecies=3)
end

# =============================================================================
# PART 1: Regression tests - ham! unchanged when couplingMat is nothing
# =============================================================================

@testset "ham! - unchanged when couplingMat is nothing (regression)" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=10, chemPot=1.0)
    b = BulkMonomer(; name="B", length=10, chemPot=1.0)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, b, -0.5)
    
    init!(l, [a, b])
    
    # Verify couplingMat is nothing
    @test a.couplingMat === nothing
    @test b.couplingMat === nothing
    
    # ham! should work normally
    ham!(a)
    ham!(b)
    
    @test a.ham !== nothing
    @test b.ham !== nothing
    @test isfinite(a.ham)
    @test isfinite(b.ham)
end

@testset "ham! - zero coupling has no effect" begin
    l = setup_lattice(15)
    
    a = BulkMonomer(; name="A", length=5, chemPot=1.0)
    b = BulkMonomer(; name="B", length=5, chemPot=1.0)
    c = BulkMonomer(; name="C", length=5, chemPot=1.0)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.3)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    init!(l, [a, b, c])
    
    ham!(a)
    ham1 = a.ham
    
    # Add zero coupling
    couplingPot!(a, b, c, 0.0)
    
    ham!(a)
    ham2 = a.ham
    
    # Should be same (zero coupling contributes nothing)
    @test isapprox(ham1, ham2, atol=1e-10)
end

# =============================================================================
# PART 2: Basic coupling contribution tests with controlled setup
# =============================================================================

@testset "ham! - coupling contributes when overlap and neighbors exist" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)  # chemPot=0 to focus on interactions
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)  # No self-interaction
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)  # No two-way
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -1.0)  # Only three-way
    
    init!(l, [a, b, c])
    
    # Clear lattices and place manually
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # A at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # C at same position (overlapping with A)
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    # B as neighbor at (6,5,5)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 6, 5, 5
    b.latS[6, 5, 5] = 1
    
    ham!(a)
    
    # Expected: coupling = nnpq(a,b,c,(5,5,5)) * -1.0
    # nnpq = 1 (b is neighbor, c overlaps)
    # ham = -1.0
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

@testset "ham! - coupling is zero when no context overlap" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear and place
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 6, 5, 5  # neighbor of A
    b.latS[6, 5, 5] = 1
    
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 10, 10  # far away, no overlap with A
    c.latS[10, 10, 10] = 1
    
    ham!(a)
    
    # Coupling should be 0 because C doesn't overlap with A
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

@testset "ham! - coupling is zero when no partner neighbors" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -1.0)
    
    init!(l, [a, b, c])
    
    # Clear and place
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5  # overlapping with A
    c.latS[5, 5, 5] = 1
    
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 10, 10, 10  # far away, not neighbor
    b.latS[10, 10, 10] = 1
    
    ham!(a)
    
    # Coupling should be 0 because B is not a neighbor of A
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

# =============================================================================
# PART 3: Manual calculation verification
# =============================================================================

@testset "ham! - matches manual calculation" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=3, chemPot=1.0)
    b = BulkMonomer(; name="B", length=3, chemPot=0.0)
    c = BulkMonomer(; name="C", length=3, chemPot=0.0)
    
    intPot!(a, -1.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, -0.2)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -0.5)
    
    init!(l, [a, b, c])
    
    # Clear and place manually
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    # A beads at (5,5,5), (6,5,5), (7,5,5) - a line
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.pos[2][1], a.pos[2][2], a.pos[2][3] = 6, 5, 5
    a.pos[3][1], a.pos[3][2], a.pos[3][3] = 7, 5, 5
    a.latS[5, 5, 5] = 1
    a.latS[6, 5, 5] = 2
    a.latS[7, 5, 5] = 3
    
    # B bead at (5,6,5) - neighbor of A[1]
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 6, 5
    b.latS[5, 6, 5] = 1
    # Other B beads far away
    b.pos[2][1], b.pos[2][2], b.pos[2][3] = 15, 15, 15
    b.pos[3][1], b.pos[3][2], b.pos[3][3] = 15, 16, 15
    b.latS[15, 15, 15] = 2
    b.latS[15, 16, 15] = 3
    
    # C bead overlapping A[1] at (5,5,5)
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    # Other C beads far away
    c.pos[2][1], c.pos[2][2], c.pos[2][3] = 15, 15, 16
    c.pos[3][1], c.pos[3][2], c.pos[3][3] = 15, 16, 16
    c.latS[15, 15, 16] = 2
    c.latS[15, 16, 16] = 3
    
    ham!(a)
    
    # Manual calculation:
    # chemPot: -(1.0) * 3 = -3.0
    #
    # Two-way interactions for A:
    # Position A[1] at (5,5,5): nnp(A,A)=1, nnp(A,B)=1, nnp(A,C)=0
    # Position A[2] at (6,5,5): nnp(A,A)=2, nnp(A,B)=0, nnp(A,C)=1 (C at (5,5,5) is neighbor!)
    # Position A[3] at (7,5,5): nnp(A,A)=1, nnp(A,B)=0, nnp(A,C)=0
    #
    # A-A: total 4 neighbors, factor 0.5: -1.0 * 4 * 0.5 = -2.0
    # A-B: total 1 neighbor: -0.2 * 1 = -0.2
    # A-C: total 1 neighbor: 0.0 * 1 = 0.0
    #
    # Three-way coupling:
    # A[1] at (5,5,5): C overlaps, B neighbor -> nnpq = 1
    # A[2] at (6,5,5): no C overlap -> nnpq = 0
    # A[3] at (7,5,5): no C overlap -> nnpq = 0
    # Total nnpq = 1, coupling = -0.5 * 1 = -0.5
    #
    # Total expected: -3.0 - 2.0 - 0.2 + 0.0 - 0.5 = -5.7
    
    expected = -3.0 - 2.0 - 0.2 - 0.5
    @test isapprox(a.ham, expected, atol=1e-10)
end

# =============================================================================
# PART 4: Multiple couplings
# =============================================================================

@testset "ham! - multiple coupling terms sum correctly" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    d = BulkMonomer(; name="D", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(d, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(a, d, 0.0)
    intPot!(b, c, 0.0)
    intPot!(b, d, 0.0)
    intPot!(c, d, 0.0)
    
    couplingPot!(a, b, c, -0.3)
    couplingPot!(a, b, d, -0.4)
    
    init!(l, [a, b, c, d])
    
    # Clear and place
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    d.latS .= 0
    
    # A at (5,5,5)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    # B neighbor at (6,5,5)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 6, 5, 5
    b.latS[6, 5, 5] = 1
    
    # C overlapping A at (5,5,5)
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    # D also overlapping A at (5,5,5)
    d.pos[1][1], d.pos[1][2], d.pos[1][3] = 5, 5, 5
    d.latS[5, 5, 5] = 1
    
    ham!(a)
    
    # Expected:
    # Coupling (b,c): nnpq(a,b,c,(5,5,5)) = 1 * -0.3 = -0.3
    # Coupling (b,d): nnpq(a,b,d,(5,5,5)) = 1 * -0.4 = -0.4
    # Total: -0.7
    
    @test isapprox(a.ham, -0.7, atol=1e-10)
end

# =============================================================================
# PART 5: Symmetry check - no double counting
# =============================================================================

@testset "ham! - no double counting of symmetric (j,k) pairs" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -1.0)  # This also sets a,c,b due to symmetry
    
    init!(l, [a, b, c])
    
    # Place B as neighbor, C as overlap
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 6, 5, 5
    b.latS[6, 5, 5] = 1
    
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 5, 5, 5
    c.latS[5, 5, 5] = 1
    
    ham!(a)
    
    # If properly avoiding double-counting, coupling contributes -1.0 once
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

@testset "ham! - symmetric configuration: B overlap, C neighbor" begin
    l = setup_lattice(20)
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    c = BulkMonomer(; name="C", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(c, 0.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -1.0)
    
    init!(l, [a, b, c])
    
    # Now B overlaps, C is neighbor (opposite of previous test)
    a.latS .= 0
    b.latS .= 0
    c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 5, 5, 5  # B overlaps
    b.latS[5, 5, 5] = 1
    
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 6, 5, 5  # C is neighbor
    c.latS[6, 5, 5] = 1
    
    ham!(a)
    
    # Due to symmetry in coupling, this should also give -1.0
    # nnpq(a, b, c, pos) counts b neighbors when c overlaps = 0
    # nnpq(a, c, b, pos) counts c neighbors when b overlaps = 1
    # But we only loop over idx_j < idx_k, so which one gets counted?
    # The coupling value is the same for both orderings, and we count
    # whichever (j,k) pair has idx_j < idx_k
    
    @test isapprox(a.ham, -1.0, atol=1e-10)
end

# =============================================================================
# PART 6: PolyMonomer support
# =============================================================================

@testset "ham! - coupling with PolyMonomer context" begin
    l = Lattice(; boxlength=20, nspecies=2)
    
    poly = Polymer(; length=5, maxlength=10, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    intPot!(a, b, 0.0)
    
    couplingPot!(a, b, pm, -0.5)  # pm is context
    
    init!(l, [poly, a, b])
    
    # Find a polymer bead position
    poly_pos = Int32[poly.pos[poly.head][1], poly.pos[poly.head][2], poly.pos[poly.head][3]]
    
    # Clear bulk lattices and place
    a.latS .= 0
    b.latS .= 0
    
    # Place A overlapping with polymer bead
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = poly_pos[1], poly_pos[2], poly_pos[3]
    a.latS[poly_pos[1], poly_pos[2], poly_pos[3]] = 1
    
    # Place B as neighbor of A
    neighbor_x = mod1(poly_pos[1] + 1, l.boxlength)
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = neighbor_x, poly_pos[2], poly_pos[3]
    b.latS[neighbor_x, poly_pos[2], poly_pos[3]] = 1
    
    ham!(a)
    
    # nnpq(a, b, pm, pos) = count b neighbors when pm overlaps at pos
    # pm overlaps at poly_pos (polymer bead there), b is neighbor
    # So coupling = -0.5 * 1 = -0.5
    @test isapprox(a.ham, -0.5, atol=1e-10)
end

# =============================================================================
# PART 7: Performance - reasonable time
# =============================================================================

@testset "ham! - performance with couplings" begin
    l = setup_lattice(30)
    
    a = BulkMonomer(; name="A", length=100, chemPot=1.0)
    b = BulkMonomer(; name="B", length=100, chemPot=1.0)
    c = BulkMonomer(; name="C", length=100, chemPot=1.0)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, -0.5)
    intPot!(b, c, -0.5)
    
    couplingPot!(a, b, c, -0.3)
    
    init!(l, [a, b, c])
    
    # Time ham! - should complete reasonably fast
    t = @elapsed for _ in 1:10
        ham!(a)
    end
    
    @test t < 5.0  # Should complete 10 calls in under 5 seconds
    @test isfinite(a.ham)
end


println()
println("All ham! coupling tests completed!")
