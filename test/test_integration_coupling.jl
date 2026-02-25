# Test file for Sections 6 & 7: FastBulkMonomer exclusion and integration tests
# for three-way coupling interactions

using Test
using Droplet

println("Testing FastBulkMonomer exclusion & integration (Sections 6-7)")
println("=" ^ 65)

# =============================================================================
# SECTION 6: FastBulkMonomer exclusion
# =============================================================================

@testset "FastBulkMonomer - couplingPot! explicit error" begin
    f = FastBulkMonomer(; name="F", length=5)
    b = BulkMonomer(; name="B", length=5)
    c = BulkMonomer(; name="C", length=5)
    
    intPot!(f, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    
    @test_throws ErrorException couplingPot!(f, b, c, -0.5)
    
    # Error message should mention FastBulkMonomer
    try
        couplingPot!(f, b, c, -0.5)
    catch e
        @test occursin("FastBulkMonomer", e.msg)
    end
end

@testset "FastBulkMonomer - ham! regression with BulkMonomer coupling present" begin
    # FastBulkMonomer should work normally even when BulkMonomer species
    # in the same simulation have three-way couplings
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    f = FastBulkMonomer(; name="F", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)
    
    intPot!(a, -0.5)
    intPot!(f, -0.5)
    intPot!(c, -0.5)
    intPot!(a, f, -0.3)
    intPot!(a, c, -0.2)
    intPot!(f, c, -0.1)
    
    # A has coupling but F does not
    couplingPot!(a, f, c, -1.0)
    
    init!(l, [a, f, c])
    
    # ham! should work on FastBulkMonomer without error
    ham!(f)
    @test isfinite(f.ham)
    
    # ham! should work on BulkMonomer with coupling
    ham!(a)
    @test isfinite(a.ham)
end

@testset "FastBulkMonomer - deltaHam!/move! regression with coupling present" begin
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=15, chemPot=0.5)
    f = FastBulkMonomer(; name="F", length=15, chemPot=0.5)
    c = BulkMonomer(; name="C", length=15, chemPot=0.5)
    
    intPot!(a, -0.5)
    intPot!(f, -0.5)
    intPot!(c, -0.5)
    intPot!(a, f, -0.3)
    intPot!(a, c, -0.2)
    intPot!(f, c, -0.1)
    
    couplingPot!(a, f, c, -0.5)
    
    init!(l, [a, f, c])
    
    ham!(a)
    ham!(f)
    ham!(c)
    
    # FastBulkMonomer move! should work normally
    n_moves = 0
    for _ in 1:200
        if move!(f)
            n_moves += 1
        end
    end
    @test n_moves > 0
    
    # After moves, FastBulkMonomer deltaH should be correct (two-way only)
    ham!(f)
    hams_f = [f.ham]
    acc_f = [f.ham]
    acc = 0.0
    
    for _ in 1:500
        if move!(f)
            acc += getDeltaHam(f)
            push!(acc_f, acc_f[1] + acc)
            ham!(f)
            push!(hams_f, f.ham)
        end
    end
    
    @test length(hams_f) > 1
    for i in 1:length(hams_f)
        @test isapprox(hams_f[i], acc_f[i], atol=1e-10)
    end
end

@testset "FastBulkMonomer - has no couplingMat field" begin
    f = FastBulkMonomer(; name="F", length=5)
    @test !hasfield(typeof(f), :couplingMat)
    @test !hasfield(typeof(f), :couplingIdx)
    @test !hasfield(typeof(f), :couplingSpecies)
end

# =============================================================================
# SECTION 7: Integration Tests
# =============================================================================

# ---- Test 1: Basic Three-Species BulkMonomer Test ----

@testset "Integration - three BulkMonomer species with coupling" begin
    l = Lattice(; boxlength=20, nspecies=3)
    
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
    
    # Manually place: A and C overlap at (10,5,5), B is neighbor at (11,5,5)
    a.latS .= 0; b.latS .= 0; c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 10, 5, 5
    a.latS[10, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 11, 5, 5
    b.latS[11, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1
    
    # ham! hand calculation:
    #   chemPot = 0, intPot = 0 for all
    #   coupling: nnpq(a,b,c,(10,5,5)) = c overlaps at (10,5,5)? yes → b neighbors? (11,5,5) is neighbor → 1
    #   coupling: nnpq(a,c,b,(10,5,5)) = b overlaps at (10,5,5)? b.latS[10,5,5]=0 → 0
    #   Total coupling = -1.0 * (1 + 0) = -1.0
    ham!(a)
    @test isapprox(a.ham, -1.0, atol=1e-10)
    
    # Now move A from (10,5,5) to (5,5,5): should lose coupling energy
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5
    
    deltaHam!(a)
    @test isapprox(a.deltaH, 1.0, atol=1e-10)  # +1.0 because lost -1.0 coupling
    
    makeMove!(a)
    ham!(a)
    @test isapprox(a.ham, 0.0, atol=1e-10)
end

# ---- Test 2: BulkMonomer + PolyMonomer coupling ----

@testset "Integration - BulkMonomer with PolyMonomer context" begin
    l = Lattice(; boxlength=20, nspecies=2)
    
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm = poly.monomers[1]
    
    a = BulkMonomer(; name="A", length=5, chemPot=0.0)
    b = BulkMonomer(; name="B", length=5, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, pm, 0.0)
    intPot!(b, pm, 0.0)
    intPot!(a, b, 0.0)
    
    couplingPot!(a, b, pm, -0.5)
    
    init!(l, [poly, a, b])
    
    # Coupling should only activate at polymer sites
    # Test consistency: ham! before + deltaH == ham! after
    n_consistent = 0
    n_trials = 100
    
    for _ in 1:n_trials
        ham!(a)
        ham_before = a.ham
        
        a.moveIdx = rand(Int32(1):Int32(a.length))
        identifyBond!(a)
        proposeMove!(a)
        
        if a.latS[a.newPos[1], a.newPos[2], a.newPos[3]] > 0
            n_consistent += 1
            continue
        end
        
        deltaHam!(a)
        delta = a.deltaH
        
        makeMove!(a)
        ham!(a)
        ham_after = a.ham
        
        if isapprox(ham_before + delta, ham_after, atol=1e-10)
            n_consistent += 1
        end
    end
    
    @test n_consistent == n_trials
end

# ---- Test 3: Heterogeneous Polymer - coupling distinguishes PolyMonomer types ----

@testset "Integration - heterogeneous polymer coupling discrimination" begin
    l = Lattice(; boxlength=20, nspecies=3)
    
    poly = Polymer(; length=10, maxlength=20, nMonomerTypes=2, pattern=:block)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
    pm1 = poly.monomers[1]  # base type
    pm2 = poly.monomers[2]  # pattern type
    
    a = BulkMonomer(; name="A", length=1, chemPot=0.0)
    b = BulkMonomer(; name="B", length=1, chemPot=0.0)
    
    intPot!(a, 0.0)
    intPot!(b, 0.0)
    intPot!(a, pm1, 0.0)
    intPot!(a, pm2, 0.0)
    intPot!(b, pm1, 0.0)
    intPot!(b, pm2, 0.0)
    intPot!(a, b, 0.0)
    
    # Coupling only with pm1, NOT pm2
    couplingPot!(a, b, pm1, -1.0)
    
    init!(l, [poly, a, b])
    
    # Find a site occupied by pm1 and a site occupied by pm2
    pm1_site = nothing
    pm2_site = nothing
    for i in 1:poly.length
        pos = poly.pos[i]
        intIdx = pos[6]  # stored intIdx
        if intIdx == pm1.intIdx && pm1_site === nothing
            pm1_site = (pos[1], pos[2], pos[3])
        elseif intIdx == pm2.intIdx && pm2_site === nothing
            pm2_site = (pos[1], pos[2], pos[3])
        end
        if pm1_site !== nothing && pm2_site !== nothing
            break
        end
    end
    
    @test pm1_site !== nothing
    @test pm2_site !== nothing
    
    # Place B as neighbor of pm1_site
    # Find a neighbor of pm1_site
    bx = mod(pm1_site[1], a.boxlength) + 1  # neighbor in x
    
    a.latS .= 0; b.latS .= 0
    
    # Place A at pm1_site (overlapping with pm1)
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm1_site[1], pm1_site[2], pm1_site[3]
    a.latS[pm1_site[1], pm1_site[2], pm1_site[3]] = 1
    
    # Place B as neighbor
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = bx, pm1_site[2], pm1_site[3]
    b.latS[bx, pm1_site[2], pm1_site[3]] = 1
    
    ham!(a)
    ham_at_pm1 = a.ham
    
    # Now place A at pm2_site instead (overlapping with pm2, not pm1)
    a.latS .= 0; b.latS .= 0
    
    bx2 = mod(pm2_site[1], a.boxlength) + 1
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = pm2_site[1], pm2_site[2], pm2_site[3]
    a.latS[pm2_site[1], pm2_site[2], pm2_site[3]] = 1
    
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = bx2, pm2_site[2], pm2_site[3]
    b.latS[bx2, pm2_site[2], pm2_site[3]] = 1
    
    ham!(a)
    ham_at_pm2 = a.ham
    
    # At pm1 site: coupling is active (-1.0 * nnpq contribution)
    # At pm2 site: coupling with pm1 should NOT be active
    # The coupling with pm1 should contribute energy at the pm1 site but not pm2
    @test ham_at_pm1 < ham_at_pm2  # pm1 site has coupling energy (negative), pm2 has 0
    @test isapprox(ham_at_pm2, 0.0, atol=1e-10)  # no coupling at pm2 site
end

# ---- Test 4: MC Consistency Test ----

@testset "Integration - MC consistency (self-moves, accumulated deltaH matches ham!)" begin
    # Coupling is "owned by self" — getDeltaHP does NOT include coupling changes
    # from partner moves. So we test self-moves only for each species, where
    # getDeltaHam correctly includes the coupling contribution.
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=30, chemPot=0.8)
    b = BulkMonomer(; name="B", length=30, chemPot=0.8)
    c = BulkMonomer(; name="C", length=30, chemPot=0.8)
    
    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)
    
    couplingPot!(a, b, c, -0.4)
    
    init!(l, [a, b, c])
    
    ham!(a)
    ham!(b)
    ham!(c)
    
    # Test each coupling owner with self-moves (same pattern as test_deltaHam_coupling.jl Part 8)
    for owner in [a, b, c]
        # Recalculate ham for all species to get clean state
        ham!(a); ham!(b); ham!(c)
        
        hams = [getHam(owner)]
        accDelta = [getHam(owner)]
        acc = 0.0
        
        niter = 500
        
        for _ in 1:niter
            if move!(owner)
                acc += getDeltaHam(owner)
                push!(accDelta, accDelta[1] + acc)
                
                ham!(owner)
                push!(hams, getHam(owner))
            end
        end
        
        @test length(hams) > 1
        for i in 1:length(hams)
            @test isapprox(hams[i], accDelta[i], atol=1e-10)
        end
    end
end

@testset "Integration - MC consistency (all species move, single-step)" begin
    # When all species move, coupling creates asymmetric deltaHP tracking.
    # Test single-step consistency: ham_before + deltaH == ham_after for each move.
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=30, chemPot=0.8)
    b = BulkMonomer(; name="B", length=30, chemPot=0.8)
    c = BulkMonomer(; name="C", length=30, chemPot=0.8)
    
    intPot!(a, -0.5)
    intPot!(b, -0.5)
    intPot!(c, -0.5)
    intPot!(a, b, -0.3)
    intPot!(a, c, -0.2)
    intPot!(b, c, -0.1)
    
    couplingPot!(a, b, c, -0.4)
    
    init!(l, [a, b, c])
    
    ham!(a); ham!(b); ham!(c)
    
    n_consistent = 0
    n_total = 0
    niter = 500
    
    for _ in 1:niter
        for (s, idx) in a.intSpecies
            ham!(s)
            ham_before = s.ham
            
            s.moveIdx = rand(Int32(1):Int32(s.length))
            identifyBond!(s)
            proposeMove!(s)
            
            if s.latS[s.newPos[1], s.newPos[2], s.newPos[3]] > 0
                n_total += 1
                n_consistent += 1
                continue
            end
            
            deltaHam!(s)
            delta = s.deltaH
            
            makeMove!(s)
            ham!(s)
            ham_after = s.ham
            
            n_total += 1
            if isapprox(ham_before + delta, ham_after, atol=1e-10)
                n_consistent += 1
            end
            
            # Recalculate all to keep state clean
            ham!(a); ham!(b); ham!(c)
        end
    end
    
    @test n_total > 0
    @test n_consistent == n_total
end

# ---- Test 5: Detailed Balance Test ----

@testset "Integration - detailed balance (forward/reverse energy symmetry)" begin
    l = Lattice(; boxlength=20, nspecies=3)
    
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
    
    # Place: A at (5,5,5), B at (11,5,5), C at (10,5,5)
    a.latS .= 0; b.latS .= 0; c.latS .= 0
    
    a.pos[1][1], a.pos[1][2], a.pos[1][3] = 5, 5, 5
    a.latS[5, 5, 5] = 1
    b.pos[1][1], b.pos[1][2], b.pos[1][3] = 11, 5, 5
    b.latS[11, 5, 5] = 1
    c.pos[1][1], c.pos[1][2], c.pos[1][3] = 10, 5, 5
    c.latS[10, 5, 5] = 1
    
    ham!(a)
    ham_initial = a.ham  # should be 0
    
    # Forward move: A from (5,5,5) to (10,5,5)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 10, 5, 5
    
    deltaHam!(a)
    delta_forward = a.deltaH
    
    # Execute forward move
    makeMove!(a)
    ham!(a)
    ham_final = a.ham  # should be -1.0
    
    @test isapprox(ham_initial + delta_forward, ham_final, atol=1e-10)
    
    # Reverse move: A from (10,5,5) back to (5,5,5)
    a.moveIdx = Int32(1)
    identifyBond!(a)
    a.newPos[1], a.newPos[2], a.newPos[3] = 5, 5, 5
    
    deltaHam!(a)
    delta_reverse = a.deltaH
    
    # Forward and reverse deltas should be equal and opposite
    @test isapprox(delta_forward + delta_reverse, 0.0, atol=1e-10)
    
    # Execute reverse move
    makeMove!(a)
    ham!(a)
    ham_restored = a.ham
    
    # Should be back to initial energy
    @test isapprox(ham_restored, ham_initial, atol=1e-10)
end

# ---- Test 6: Multiple coupling owners ----

@testset "Integration - multiple BulkMonomers each owning couplings" begin
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    c = BulkMonomer(; name="C", length=20, chemPot=0.5)
    
    intPot!(a, -0.3)
    intPot!(b, -0.3)
    intPot!(c, -0.3)
    intPot!(a, b, -0.2)
    intPot!(a, c, -0.1)
    intPot!(b, c, -0.1)
    
    # Both A and B own couplings
    couplingPot!(a, b, c, -0.5)
    couplingPot!(b, a, c, -0.3)
    
    init!(l, [a, b, c])
    
    # Both should have coupling matrices
    @test a.couplingMat !== nothing
    @test b.couplingMat !== nothing
    @test c.couplingMat === nothing  # c doesn't own any coupling
    
    # Consistency test for both A and B
    for owner in [a, b]
        n_consistent = 0
        n_trials = 100
        
        for _ in 1:n_trials
            ham!(owner)
            ham_before = owner.ham
            
            owner.moveIdx = rand(Int32(1):Int32(owner.length))
            identifyBond!(owner)
            proposeMove!(owner)
            
            if owner.latS[owner.newPos[1], owner.newPos[2], owner.newPos[3]] > 0
                n_consistent += 1
                continue
            end
            
            deltaHam!(owner)
            delta = owner.deltaH
            
            makeMove!(owner)
            ham!(owner)
            ham_after = owner.ham
            
            if isapprox(ham_before + delta, ham_after, atol=1e-10)
                n_consistent += 1
            end
        end
        
        @test n_consistent == n_trials
    end
end

# ---- Test 7: Mixed FastBulkMonomer + BulkMonomer simulation ----

@testset "Integration - mixed FastBulkMonomer and BulkMonomer simulation" begin
    l = Lattice(; boxlength=15, nspecies=3)
    
    a = BulkMonomer(; name="A", length=20, chemPot=0.8)
    f = FastBulkMonomer(; name="F", length=20, chemPot=0.8)
    c = BulkMonomer(; name="C", length=20, chemPot=0.8)
    
    intPot!(a, -0.5)
    intPot!(f, -0.5)
    intPot!(c, -0.5)
    intPot!(a, f, -0.3)
    intPot!(a, c, -0.2)
    intPot!(f, c, -0.1)
    
    couplingPot!(a, f, c, -0.4)
    
    init!(l, [a, f, c])
    
    # A's move consistency
    ham!(a); ham!(f); ham!(c)
    
    a_hams = [getHam(a)]
    a_accDelta = [getHam(a)]
    a_acc = 0.0
    a.deltaH = 0.0; a.deltaHP = 0.0
    
    niter = 500
    for _ in 1:niter
        if move!(a)
            a_acc += getDeltaHam(a)
            push!(a_accDelta, a_accDelta[1] + a_acc)
            ham!(a)
            push!(a_hams, getHam(a))
        end
    end
    
    @test length(a_hams) > 1
    for i in 1:length(a_hams)
        @test isapprox(a_hams[i], a_accDelta[i], atol=1e-10)
    end
    
    # F's move consistency (should be pure two-way, no coupling)
    ham!(f)
    f_hams = [getHam(f)]
    f_accDelta = [getHam(f)]
    f_acc = 0.0
    f.deltaH = 0.0; f.deltaHP = 0.0
    
    for _ in 1:niter
        if move!(f)
            f_acc += getDeltaHam(f)
            push!(f_accDelta, f_accDelta[1] + f_acc)
            ham!(f)
            push!(f_hams, getHam(f))
        end
    end
    
    @test length(f_hams) > 1
    for i in 1:length(f_hams)
        @test isapprox(f_hams[i], f_accDelta[i], atol=1e-10)
    end
end

# ---- Test 8: BulkMonomer + Polymer full simulation ----

@testset "Integration - BulkMonomer coupling with Polymer in simulation" begin
    l = Lattice(; boxlength=15, nspecies=3)
    
    poly = Polymer(; length=8, maxlength=16, nMonomerTypes=1)
    initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=-0.2)
    pm = poly.monomers[1]
    
    a = BulkMonomer(; name="A", length=20, chemPot=0.5)
    b = BulkMonomer(; name="B", length=20, chemPot=0.5)
    
    intPot!(a, -0.3)
    intPot!(b, -0.3)
    intPot!(a, pm, -0.1)
    intPot!(b, pm, -0.1)
    intPot!(a, b, -0.2)
    
    couplingPot!(a, b, pm, -0.5)
    
    init!(l, [poly, a, b])
    
    # A's self-move consistency with polymer present
    ham!(a); ham!(b)
    ham!(pm)
    
    a_hams = [getHam(a)]
    a_accDelta = [getHam(a)]
    a_acc = 0.0
    
    niter = 500
    for _ in 1:niter
        if move!(a)
            a_acc += getDeltaHam(a)
            push!(a_accDelta, a_accDelta[1] + a_acc)
            ham!(a)
            push!(a_hams, getHam(a))
        end
    end
    
    @test length(a_hams) > 1
    for i in 1:length(a_hams)
        @test isapprox(a_hams[i], a_accDelta[i], atol=1e-10)
    end
end

# ---- Test 9: Performance ----

@testset "Integration - performance acceptable" begin
    l = Lattice(; boxlength=30, nspecies=3)
    
    a = BulkMonomer(; name="A", length=200, chemPot=1.0)
    b = BulkMonomer(; name="B", length=200, chemPot=1.0)
    c = BulkMonomer(; name="C", length=200, chemPot=1.0)
    
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, -0.5)
    intPot!(b, c, -0.5)
    
    couplingPot!(a, b, c, -0.3)
    
    init!(l, [a, b, c])
    
    # Warmup
    ham!(a)
    for _ in 1:10
        move!(a)
    end
    
    # Time 1000 move attempts (includes propose + deltaHam + overlap check + makeMove)
    t = @elapsed for _ in 1:1000
        move!(a)
    end
    
    @test t < 10.0  # Should be well under 10 seconds
    @test isfinite(a.ham)
    
    println("  Performance: 1000 move! calls in $(round(t, digits=3))s")
end

println()
println("All Section 6-7 tests completed!")
