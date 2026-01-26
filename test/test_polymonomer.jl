# MARK: -PolyMonomer init
@testset "Basic PolyMonomer initialization" begin
    # Create a simple polymer with 2 monomer types
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    # Initialize monomers with default potentials
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    # Check that monomers were created
    @test poly.monomers !== nothing
    @test length(poly.monomers) == 2
    # Note: monomer types are now stored in pos[key][6] after init!() is called,
    # and patternTails tracks pattern regions instead of monomerIndex.
    # Before init!(), patternTails is nothing since types aren't assigned to pos yet.
    
    # Check that each monomer has the poly reference
    @test poly.monomers[1].poly === poly
    @test poly.monomers[2].poly === poly
    
    # Check that monomers inherit chemical and interaction potentials
    @test poly.monomers[1].chemPot == 0.0
    @test poly.monomers[2].chemPot == 0.0
    @test poly.monomers[1].intPot == 0.0
    @test poly.monomers[2].intPot == 0.0
    
    # Check that interaction structures are initialized
    @test poly.monomers[1].intSpecies !== nothing
    @test poly.monomers[1].intMat !== nothing
    @test poly.monomers[1].intIdx !== nothing
    
    # Check that both monomers share the same interaction matrix
    @test poly.monomers[1].intMat === poly.monomers[2].intMat
    
    # Now initialize with lattice to assign types
    l = Lattice(; boxlength = 8, nspecies = 1)
    init!(l, [poly])
    
    # After init!(), patternTails should be populated for polymers with patterns
    @test poly.patternTails !== nothing
end

# MARK: -PolyMonomer pattern
@testset "PolyMonomer pattern assignment" begin
    # Test uniform pattern (default) - requires init!() to assign types to pos
    l1 = Lattice(; boxlength = 8, nspecies = 1)
    poly1 = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :uniform
    )
    initMonomers!(poly1, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l1, [poly1])
    # All beads should have type 1 for uniform pattern
    types1 = [poly1.pos[k][6] for k in keys(poly1.pos)]
    @test all(types1 .== 1)
    
    # Test block pattern - types 2+ must be contiguous blocks
    l3 = Lattice(; boxlength = 8, nspecies = 1)
    poly3 = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    initMonomers!(poly3, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l3, [poly3])
    # Block pattern should have first half type 1, second half type 2
    types3 = [poly3.pos[k][6] for k in keys(poly3.pos)]
    @test count(types3 .== 1) == 2  # Half of 4 beads
    @test count(types3 .== 2) == 2  # Other half
    
    # Test custom pattern with contiguous blocks
    l4 = Lattice(; boxlength = 8, nspecies = 1)
    poly4 = Polymer(;
        length = 8,
        maxlength = 8,
        nMonomerTypes = 3,
        pattern = :custom,
        customPattern = [1, 1, 2, 2, 3, 3, 1, 1]  # Contiguous blocks
    )
    initMonomers!(poly4, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l4, [poly4])
    
    # Verify types are assigned correctly by walking the chain
    types4 = Int32[]
    key = poly4.head
    while key != Int32(0)
        push!(types4, poly4.pos[key][6])
        key = poly4.pos[key][5]  # next
    end
    @test types4 == Int32[1, 1, 2, 2, 3, 3, 1, 1]
end

# MARK: -PolyMonomer intPot
@testset "Setting interaction potentials between monomers" begin
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 3
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    activator = poly.monomers[3]
    
    dna.name = "dna"
    promoter.name = "promoter"
    activator.name = "activator"
    
    # Set self-interaction
    intPot!(dna, 0.5)
    @test dna.intPot == 0.5
    @test dna.intMat[dna.intIdx, dna.intIdx] == 0.5
    
    # Set cross-interaction
    intPot!(promoter, activator, 2.0)
    @test promoter.intMat[promoter.intIdx, activator.intIdx] == 2.0
    @test activator.intMat[activator.intIdx, promoter.intIdx] == 2.0
    
    # Check that all monomers still share the same matrix
    @test dna.intMat === promoter.intMat
    @test promoter.intMat === activator.intMat
end

# MARK: -PolyMono FastBulk intPot
@testset "Interaction between FastBulkMonomer and PolyMonomer" begin
    l = Lattice(;
        boxlength = 4,
        nspecies = 3
    )
    
    tf = FastBulkMonomer(;
        length = 12,
        initShape = "uniform"
    )
    
    poly = Polymer(;
        length = l.boxlength,
        maxlength = 8,
        nMonomerTypes = 3
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    
    dna.name = "dna"
    promoter.name = "promoter"
    
    # Set interaction potentials between bulk monomer and poly monomers
    intPot!(tf, dna, 0.3)
    intPot!(tf, promoter, 2.0)
    
    # Check that interactions were set correctly
    @test getIntPot(tf, dna) == 0.3
    @test getIntPot(dna, tf) == 0.3
    @test getIntPot(tf, promoter) == 2.0
    @test getIntPot(promoter, tf) == 2.0
    
    # Check that all share the same interaction matrix
    @test tf.intMat === dna.intMat
    @test dna.intMat === promoter.intMat
end

# MARK: -PolyMonomer mixed
@testset "Full initialization with mixed species" begin
    l = Lattice(;
        boxlength = 8,
        nspecies = 3
    )
    
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block  # Changed from :alternating to :block for contiguous blocks
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    
    dna.name = "dna"
    promoter.name = "promoter"
    
    intPot!(tf, 1.0)
    intPot!(tf, dna, 0.3)
    intPot!(tf, promoter, 2.0)
    
    # Initialize the species array with mixed types
    species = [tf, dna, promoter]
    
    # This should work with the new init! implementation
    init!(l, species)
    
    # Verify lattice was initialized
    @test l.lat !== nothing
    @test l.nspecies > 0
    
    # Verify polymer has lattice view
    @test poly.latS !== nothing
    @test poly.boxlength == 8
    
    # Verify bulk monomer has lattice view
    @test tf.latS !== nothing
    @test tf.boxlength == 8
end

# MARK: -PolyMonomer helpers
@testset "Position and lattice delegation helpers" begin
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    l = Lattice(;
        boxlength = 8,
        nspecies = 1
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    monomer = poly.monomers[1]
    
    # Test delegation helpers
    @test getLatS(monomer) === poly.latS
    @test getTempVec(monomer) === poly.tempVec
    @test getBoxlength(monomer) == poly.boxlength
    
    # Test with regular species for comparison
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    l2 = Lattice(;
        boxlength = 8,
        nspecies = 1
    )
    
    init!(l2, [tf])
    
    @test getLatS(tf) === tf.latS
    @test getTempVec(tf) === tf.tempVec
    @test getBoxlength(tf) == tf.boxlength
end

# MARK: -PolyMonomer nn/nnp
@testset "nn and nnp with PolyMonomer" begin
    l = Lattice(;
        boxlength = 8,
        nspecies = 2
    )
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly, tf])
    
    monomer = poly.monomers[1]
    
    # Test nn (same-species neighbors)
    pos = poly.pos[1][1:3]
    neighbors = nn(monomer, pos)
    @test neighbors >= 0
    @test neighbors <= 6
    
    # Test nnp (cross-species neighbors)
    neighbors_cross = nnp(monomer, tf, pos)
    @test neighbors_cross >= 0
    @test neighbors_cross <= 6
end