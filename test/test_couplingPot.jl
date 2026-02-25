# Test file for Section 3: couplingPot! function
# Tests the three-way coupling initialization interface

using Test
using Droplet

println("Testing couplingPot! (Section 3)")
println("=" ^ 50)

# =============================================================================
# PART 1: Basic functionality
# =============================================================================

@testset "couplingPot! - basic setting and retrieval" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    # Set up two-way interactions first
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    # Set three-way coupling
    couplingPot!(a, b, c, -0.5)
    
    # Verify coupling was set
    @test getCouplingPot(a, b, c) == -0.5
end

@testset "couplingPot! - symmetry in j,k enforced" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    couplingPot!(a, b, c, -0.5)
    
    # Both orderings should return the same value
    @test getCouplingPot(a, b, c) == getCouplingPot(a, c, b)
    @test getCouplingPot(a, b, c) == -0.5
    @test getCouplingPot(a, c, b) == -0.5
end

# =============================================================================
# PART 2: Validation - distinct species required
# =============================================================================

@testset "couplingPot! - error when i === j" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    
    @test_throws ErrorException couplingPot!(a, a, b, -0.5)
end

@testset "couplingPot! - error when j === k" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    
    @test_throws ErrorException couplingPot!(a, b, b, -0.5)
end

@testset "couplingPot! - error when i === k" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    
    @test_throws ErrorException couplingPot!(a, b, a, -0.5)
end

# =============================================================================
# PART 3: Auto-initialization of structures
# =============================================================================

@testset "couplingPot! - auto-initializes intSpecies" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    # Don't call intPot! first - couplingPot! should auto-initialize
    @test a.intSpecies === nothing
    @test b.intSpecies === nothing
    @test c.intSpecies === nothing
    
    couplingPot!(a, b, c, -0.5)
    
    # Now all should have intSpecies
    @test a.intSpecies !== nothing
    @test b.intSpecies !== nothing
    @test c.intSpecies !== nothing
    @test a.intIdx !== nothing
    @test b.intIdx !== nothing
    @test c.intIdx !== nothing
end

@testset "couplingPot! - initializes couplingIdx and couplingSpecies" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    @test a.couplingIdx === nothing
    @test a.couplingSpecies === nothing
    @test a.couplingMat === nothing
    
    couplingPot!(a, b, c, -0.5)
    
    @test a.couplingIdx !== nothing
    @test a.couplingSpecies !== nothing
    @test a.couplingMat !== nothing
    @test a.couplingIdx == 1  # First BulkMonomer in coupling list
end

# =============================================================================
# PART 4: Multiple couplings and species
# =============================================================================

@testset "couplingPot! - multiple (j,k) pairs work correctly" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    d = BulkMonomer(; name="D", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(d, -1.0)
    
    couplingPot!(a, b, c, -0.5)
    couplingPot!(a, b, d, -0.3)
    couplingPot!(a, c, d, -0.2)
    
    @test getCouplingPot(a, b, c) == -0.5
    @test getCouplingPot(a, b, d) == -0.3
    @test getCouplingPot(a, c, d) == -0.2
    
    # Symmetry still holds
    @test getCouplingPot(a, c, b) == -0.5
    @test getCouplingPot(a, d, b) == -0.3
    @test getCouplingPot(a, d, c) == -0.2
end

@testset "couplingPot! - multiple BulkMonomers can have independent couplings" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    # A has coupling with B,C
    couplingPot!(a, b, c, -0.5)
    
    # B has different coupling with A,C
    couplingPot!(b, a, c, -0.8)
    
    @test getCouplingPot(a, b, c) == -0.5
    @test getCouplingPot(b, a, c) == -0.8
    
    # They should be independent
    @test getCouplingPot(a, b, c) != getCouplingPot(b, a, c)
end

@testset "couplingPot! - adding new species expands matrices correctly" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    couplingPot!(a, b, c, -0.5)
    
    # Now add a new species D
    d = BulkMonomer(; name="D", length=5)
    intPot!(d, -1.0)
    intPot!(a, d, 0.0)  # This should expand intMat
    
    # Original coupling should still work
    @test getCouplingPot(a, b, c) == -0.5
    
    # Can add new coupling involving D
    couplingPot!(a, b, d, -0.3)
    
    @test getCouplingPot(a, b, d) == -0.3
    @test getCouplingPot(a, b, c) == -0.5  # Original still preserved
end

@testset "couplingPot! - coupling values survive matrix sync" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    # Set up A's couplings
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    couplingPot!(a, b, c, -0.5)
    
    # Create independent species and link them
    d = BulkMonomer(; name="D", length=5)
    e = BulkMonomer(; name="E", length=5)
    intPot!(d, -1.0)
    intPot!(e, -1.0)
    intPot!(d, e, 0.0)
    
    # Now link the networks by setting interaction between a and d
    intPot!(a, d, 0.0)
    
    # Original coupling should be preserved
    @test getCouplingPot(a, b, c) == -0.5
end

# =============================================================================
# PART 5: PolyMonomer support
# =============================================================================

@testset "couplingPot! - PolyMonomer as partner" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    pm = poly.monomers[1]
    
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    
    # PolyMonomer as partner (j)
    couplingPot!(a, pm, b, -0.5)
    
    @test getCouplingPot(a, pm, b) == -0.5
    @test getCouplingPot(a, b, pm) == -0.5  # Symmetry
end

@testset "couplingPot! - PolyMonomer as context" begin
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.5)
    pm = poly.monomers[1]
    
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    
    # PolyMonomer as context (k)
    couplingPot!(a, b, pm, -0.5)
    
    @test getCouplingPot(a, b, pm) == -0.5
    @test getCouplingPot(a, pm, b) == -0.5  # Symmetry
end

# =============================================================================
# PART 6: Zero coupling and edge cases
# =============================================================================

@testset "couplingPot! - zero coupling can be set" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    couplingPot!(a, b, c, 0.0)
    
    @test getCouplingPot(a, b, c) == 0.0
end

@testset "couplingPot! - overwriting existing coupling" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    couplingPot!(a, b, c, -0.5)
    @test getCouplingPot(a, b, c) == -0.5
    
    # Overwrite
    couplingPot!(a, b, c, -0.8)
    @test getCouplingPot(a, b, c) == -0.8
end

@testset "getCouplingPot - returns nothing when not initialized" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    @test getCouplingPot(a, b, c) === nothing
end

@testset "getCouplingPot - returns nothing for unknown species" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    d = BulkMonomer(; name="D", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    couplingPot!(a, b, c, -0.5)
    
    # d is not in the interaction network
    @test getCouplingPot(a, b, d) === nothing
    @test getCouplingPot(a, d, c) === nothing
end

# =============================================================================
# PART 7: FastBulkMonomer exclusion
# =============================================================================

@testset "couplingPot! - error for FastBulkMonomer as self" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(f, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    # FastBulkMonomer as self should fail (MethodError because couplingPot! is only defined for BulkMonomer)
    @test_throws MethodError couplingPot!(f, b, c, -0.5)
end

@testset "couplingPot! - FastBulkMonomer allowed as partner or context" begin
    a = BulkMonomer(; name="A", length=5)
    f = FastBulkMonomer(; name="F", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(f, -1.0)
    intPot!(c, -1.0)
    
    # FastBulkMonomer as partner - should work
    couplingPot!(a, f, c, -0.5)
    @test getCouplingPot(a, f, c) == -0.5
end

# =============================================================================
# PART 8: Matrix dimensions
# =============================================================================

@testset "couplingPot! - couplingMat has correct dimensions" begin
    a = BulkMonomer(; name="A", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, 0.0)
    intPot!(a, c, 0.0)
    intPot!(b, c, 0.0)
    
    couplingPot!(a, b, c, -0.5)
    
    # couplingMat should be m × n × n where:
    # m = number of BulkMonomers in couplingSpecies (1, just a)
    # n = number of species in intMat (3: a, b, c)
    @test size(a.couplingMat, 1) == 1  # Only 'a' has couplings
    @test size(a.couplingMat, 2) == size(a.intMat, 1)
    @test size(a.couplingMat, 3) == size(a.intMat, 1)
end


println()
println("All couplingPot! tests completed!")
