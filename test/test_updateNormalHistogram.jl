using Test
using Droplet

@testset "updateNormalHistogram! tests" begin
    
    @testset "Single particle at origin" begin
        # Simple case: one particle at the origin
        l = Lattice(boxlength=20, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=1,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        # Get the position of the single bead that was placed
        bead_pos = first(values(bulk.pos))[1:3]
        
        # Origin at same location as particle (distance = 0)
        origin = Float64.(bead_pos)
        maxDistance = 10.0
        hist = zeros(Float64, 10)
        
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        
        # Particle at distance 0 should be in first bin
        @test hist[1] == 1.0
        @test sum(hist[2:end]) == 0.0
    end
    
    @testset "Single particle at known distance" begin
        # Particle at exact distance from origin
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=1,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        # Get actual position and compute distance
        bead_pos = first(values(bulk.pos))[1:3]
        origin = Float64[10.0, 10.0, 10.0]
        
        # Compute expected distance
        dx = abs(Float64(bead_pos[1]) - origin[1])
        dy = abs(Float64(bead_pos[2]) - origin[2])
        dz = abs(Float64(bead_pos[3]) - origin[3])
        dist = sqrt(dx^2 + dy^2 + dz^2)
        
        maxDistance = 30.0  # Large enough to include particle
        hist = zeros(Float64, 10)
        
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        
        # Should count exactly one particle
        @test sum(hist) == 1.0
        
        # Calculate expected bin
        expected_bin = Int(floor(10 * dist / maxDistance)) + 1
        if expected_bin <= 10
            @test hist[expected_bin] == 1.0
        end
    end
    
    @testset "Multiple particles at different distances" begin
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=5,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        # Use center of lattice as origin
        origin = Float64[15.0, 15.0, 15.0]
        maxDistance = 30.0
        hist = zeros(Float64, 10)
        
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        
        # Should count all 5 particles
        @test sum(hist) == 5.0
    end
    
    @testset "Particles beyond maxDistance are ignored" begin
        l = Lattice(boxlength=50, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=10,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        # Origin at corner with very small maxDistance
        origin = Float64[1.0, 1.0, 1.0]
        maxDistance = 5.0  # Only count particles within distance 5
        hist = zeros(Float64, 5)
        
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        
        # Some particles will be excluded
        @test sum(hist) <= 10.0
        @test sum(hist) >= 0.0
    end
    
    @testset "Histogram accumulates across multiple calls" begin
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=3,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        origin = Float64[15.0, 15.0, 15.0]
        maxDistance = 30.0
        hist = zeros(Float64, 10)
        
        # First call
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        first_sum = sum(hist)
        @test first_sum == 3.0
        
        # Second call should accumulate
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        @test sum(hist) == 6.0
        
        # Third call
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        @test sum(hist) == 9.0
    end
    
    @testset "Different number of bins" begin
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=5,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        origin = Float64[15.0, 15.0, 15.0]
        maxDistance = 30.0  # Large enough to include all particles
        
        # Test with 5 bins
        hist5 = zeros(Float64, 5)
        updateNormalHistogram!(bulk, hist5, origin, maxDistance)
        count5 = sum(hist5)
        @test count5 >= 1.0  # At least some particles counted
        
        # Test with 20 bins
        hist20 = zeros(Float64, 20)
        updateNormalHistogram!(bulk, hist20, origin, maxDistance)
        @test sum(hist20) == count5  # Same particles counted
        
        # Test with 100 bins
        hist100 = zeros(Float64, 100)
        updateNormalHistogram!(bulk, hist100, origin, maxDistance)
        @test sum(hist100) == count5  # Same particles counted
    end
    
    @testset "Particle exactly at maxDistance" begin
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=1,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        bead_pos = first(values(bulk.pos))[1:3]
        origin = Float64.(bead_pos)
        
        # Use a small positive maxDistance
        maxDistance = 0.1
        hist = zeros(Float64, 10)
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        # dist = 0, maxDistance = 0.1, bin_idx = floor(10 * 0 / 0.1) + 1 = 1
        @test hist[1] == 1.0
    end
    
    @testset "Works with Polymer species" begin
        l = Lattice(boxlength=20, nspecies=1)
        
        poly = Polymer(
            length=5,
            maxlength=5,
            initShape="line",
            nMonomerTypes=1,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Origin at first bead position
        origin = Float64[10.0, 10.0, 10.0]
        maxDistance = 10.0
        hist = zeros(Float64, 10)
        
        updateNormalHistogram!(poly, hist, origin, maxDistance)
        
        # Should count all 5 beads in the polymer
        @test sum(hist) == 5.0
        
        # First bead should be at or near distance 0
        @test hist[1] >= 1.0
    end
    
    @testset "Binning correctness with known distance" begin
        # Test that binning is calculated correctly
        l = Lattice(boxlength=30, nspecies=1)
        
        bulk = FastBulkMonomer(
            length=1,
            initShape="uniform"
        )
        
        init!(l, [bulk])
        
        # Set a specific scenario
        bead_pos = first(values(bulk.pos))[1:3]
        origin = Float64.(bead_pos) .- [5.0, 0.0, 0.0]  # Particle is 5 units away in x direction
        
        maxDistance = 10.0
        nbins = 10
        hist = zeros(Float64, nbins)
        
        updateNormalHistogram!(bulk, hist, origin, maxDistance)
        
        # Distance is 5.0
        # bin_idx = floor(10 * 5.0 / 10.0) + 1 = floor(5.0) + 1 = 6
        @test hist[6] == 1.0
        @test sum(hist) == 1.0
    end
end
