let 
    using Droplet
    using StatsBase
    using Distributions
end

## 1010011111110000100000101010111100010100001100101110000
# Saving and loading
let 
    date = "20250718"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 16,
        nspecies = 1
    )

    initLatticeRandom!(l)
    export!(l, outpath)
    save!(l, outpath)

    l0 = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    load!(l0, outpath)
    # @assert l.lat == l0.lat

end

## 1010011111110000100000101010111100010100001100101110000
# Adding species to a lattice, initializing
let 
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])
    
    for i in 1:16
        @assert s1.latS[i,8,8] == 1
        @assert l.lat[1,i,8,8] == 1
    end

    @assert s1.pos[1][4] == 0
    @assert s1.pos[5][4] == 4
    @assert s1.pos[5][5] == 6
    @assert s1.pos[16][5] == 0

    # getting information about a bond in prep for a move
    s1.moveIdx = 5

    identifyBond!(s1)

    @assert s1.ki == 5
    @assert s1.kj == 6
    @assert s1.kprev == 4
    @assert s1.knext == 7
    @assert s1.bi == [5,8,8]
    @assert s1.bj == [6,8,8]
    @assert s1.bprev == [4,8,8]
    @assert s1.bnext == [7,8,8]

    s1.moveIdx = 1

    identifyBond!(s1)

    @assert s1.ki == 1
    @assert s1.kj == 2
    @assert s1.kprev == 0

    s1.moveIdx = 15

    identifyBond!(s1)

    @assert s1.ki == 15
    @assert s1.kj == 16
    @assert s1.knext == 0

    export!(l, outpath)
end

## 1010011111110000100000101010111100010100001100101110000
# Moving bonds
let
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])

    ### ### ### Adding a bond

    s1.moveIdx = 5
    identifyBond!(s1)
    s1.moveType = 1  # bond addition
    s1.moveVec = [0,0,1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    # Original length
    origLength = s1.length
    
    makeMove!(s1)
    
    # Tests
    @assert s1.length == origLength + 2  # Length increased by 2
    @assert haskey(s1.pos, 17)  # New positions exist
    @assert haskey(s1.pos, 18)
    @assert s1.pos[5][5] == 17  # Check connections were updated
    @assert s1.pos[17][5] == 18
    @assert s1.pos[18][5] == 6
    @assert s1.pos[18][4] == 17
    @assert s1.pos[17][4] == 5
    @assert s1.pos[5][4] == 4
    @assert s1.latS[5,8,9] == 1  # New lattice sites are occupied
    @assert s1.latS[6,8,9] == 1

    ### ### ### Kinking a bond

    s1.moveIdx = 4
    identifyBond!(s1)
    s1.moveType = 2  # kinks
    s1.moveVec = [0,0,1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    # Original length
    origLength = s1.length
    
    makeMove!(s1)
    
    # Tests
    @assert s1.length == origLength  # Length unchanged
    @assert s1.pos[5][5] == 17  # Check connections remain the same
    @assert s1.pos[17][5] == 18
    @assert s1.pos[18][5] == 6
    @assert s1.pos[18][4] == 17
    @assert s1.pos[17][4] == 5
    @assert s1.pos[5][4] == 4
    @assert s1.latS[4,8,9] == 1  # New lattice site is occupied
    @assert s1.latS[5,8,8] == 0  # Old lattice site is zero

    s1.moveIdx = 17
    identifyBond!(s1)
    s1.moveType = 2  # kinks
    s1.moveVec = [0,0,-1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    # Original length
    origLength = s1.length
    
    makeMove!(s1)
    
    # Tests
    @assert s1.length == origLength  # Length unchanged
    @assert s1.pos[5][5] == 17  # Check connections remain the same
    @assert s1.pos[17][5] == 18
    @assert s1.pos[18][5] == 6
    @assert s1.pos[18][4] == 17
    @assert s1.pos[17][4] == 5
    @assert s1.pos[5][4] == 4
    @assert s1.latS[6,8,9] == 0  # New lattice site is occupied
    @assert s1.latS[6,8,8] == 1  # Old lattice site is zero

    ### ### ### Removing a bond

    s1.moveIdx = 5
    identifyBond!(s1)
    s1.moveType = 3  # remove
    s1.moveVec = [0,0,-1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    # Original length
    origLength = s1.length
    
    makeMove!(s1)
    
    # Tests
    @assert s1.length == origLength - 2  # Length reduced by 2
    @assert s1.pos[4][5] == 18  # Check connections updated
    @assert s1.pos[18][5] == 6
    @assert s1.pos[6][4] == 18
    @assert s1.pos[18][4] == 4
    @assert s1.latS[5,8,9] == 0  
    @assert s1.latS[5,8,8] == 1

    # v = view(s1.pos[18], 1:3)
    # s1.latS[v...] = 0

    export!(l, outpath)

end

## 1010011111110000100000101010111100010100001100101110000
# Testing moveTypes!
let
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])

    ### ### ### Adding a bond

    s1.moveIdx = 5
    identifyBond!(s1)
    s1.moveType = 1  # bond addition
    s1.moveVec = [0,0,1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    # Original length
    origLength = s1.length
    
    makeMove!(s1)

    s1.moveIdx = 3
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [1,1,1,1]

    s1.moveIdx = 4
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [1,1,2,1]

    s1.moveIdx = 5
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [2,2,1,1]

    s1.moveIdx = 17
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [1,1,1,3]

    s1.moveIdx = 18
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [2,2,1,1]

    
    s1.moveIdx = 18
    identifyBond!(s1)
    s1.moveType = 1  # bond addition
    s1.moveVec = [0,1,0]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    makeMove!(s1)

    s1.moveIdx = 18
    identifyBond!(s1)
    moveTypes!(s1)

    @assert s1.moveTypes == [1,2,1,2]

    export!(l, outpath)
end

## 1010011111110000100000101010111100010100001100101110000
# Testing proposeMove!
let
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])

    ### ### ### Adding a bond

    s1.moveIdx = 5
    identifyBond!(s1)
    s1.moveType = 1  # bond addition
    s1.moveVec = [0,0,1]
    s1.newPosi = s1.bi .+ s1.moveVec
    s1.newPosj = s1.bj .+ s1.moveVec
    
    makeMove!(s1);

    s1.moveIdx = 17
    identifyBond!(s1)
    moveTypes!(s1)

    @show s1.moveTypes

    proposeMove!(s1)

    @show s1.moveType

    println("")

end

## 1010011111110000100000101010111100010100001100101110000
# Testing hamiltonian functions
let
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 16,
        maxlength = 40,
        boxlength = l.boxlength,
        intPot = 1.0,
        chemPot = 1.0
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength,
        intPot = 0.0,
        chemPot = 0.0
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])
    
    # initialize hamiltonians
    initHamPoly!(s1)
    hamInitial = s1.ham

    #### iterating multiple moves
    
    nmoves = 0

    for i in 1:100
        nmoves += move!(s1)
    end

    @show nmoves

    # ### ### ### Adding a bond

    # s1.moveIdx = 5
    # identifyBond!(s1)
    # s1.moveType = 1  # bond addition
    # s1.moveVec = [0,0,1]
    # s1.newPosi = s1.bi .+ s1.moveVec
    # s1.newPosj = s1.bj .+ s1.moveVec

    # deltaHamPoly!(s1)
    # @show s1.deltaH
    
    # makeMove!(s1)

    # hamPoly!(s1)
    # @show s1.ham
    # @show hamInitial + s1.deltaH
    # @assert s1.ham == hamInitial + s1.deltaH

    # ### ### ### Removing a bond

    # s1.moveIdx = 17
    # identifyBond!(s1)
    # s1.moveType = 3  # bond addition
    # s1.moveVec = [0,0,-1]
    # s1.newPosi = s1.bi .+ s1.moveVec
    # s1.newPosj = s1.bj .+ s1.moveVec

    # deltaHamPoly!(s1)
    # @show s1.deltaH
    
    # makeMove!(s1)

    # hamPoly!(s1)
    # @show s1.ham
    # @show hamInitial + s1.deltaH
    # @assert s1.ham == hamInitial + s1.deltaH

    export!(l, outpath)
end

## 1010011111110000100000101010111100010100001100101110000
# ensures positions are periodic within the lattice's bounding box
function boxmod!(v::Vector{Int32}, w::Vector{Int32}, lb::Int32, ub::Int32)

    for i in eachindex(v)
        if v[i] > ub
            w[i] = v[i] - (ub - lb + 1)
        elseif v[i] < lb
            w[i] = v[i] + (ub - lb + 1)
        else
            w[i] = v[i]
        end

    return w
end

let 


end

## 1010011111110000100000101010111100010100001100101110000
# What is context?

let 
    filename = "ctx.json"
    open(filename, "w") do f
        JSON.print(f, CONTEXT, 4)
    end
end

# v = view(s1.pos[18], 1:3)
# s1.latS[v...] = 0