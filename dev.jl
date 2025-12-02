let 
    using Droplet
    using StatsBase
    using Distributions
    using Plots
    using Profile
    using OrderedCollections
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
    global p = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # Bulk monomers
    global b = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(p)
    initBulk!(b)
    initLattice!(l, [p, b])
    
    for i in 1:16
        @assert p.latS[i,8,8] == 1
        @assert l.lat[1,i,8,8] == 1
    end

    @assert p.pos[1][4] == 0
    @assert p.pos[5][4] == 4
    @assert p.pos[5][5] == 6
    @assert p.pos[16][5] == 0

    # getting information about a bond in prep for a move
    p.moveIdx = 5

    identifyBond!(p)

    @assert p.ki == 5
    @assert p.kj == 6
    @assert p.kprev == 4
    @assert p.knext == 7
    @assert p.bi == [5,8,8]
    @assert p.bj == [6,8,8]
    @assert p.bprev == [4,8,8]
    @assert p.bnext == [7,8,8]

    p.moveIdx = 1

    identifyBond!(p)

    @assert p.ki == 1
    @assert p.kj == 2
    @assert p.kprev == 0

    p.moveIdx = 15

    identifyBond!(p)

    @assert p.ki == 15
    @assert p.kj == 16
    @assert p.knext == 0

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
    global p = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # Bulk monomers
    global b = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(p)
    initBulk!(b)
    initLattice!(l, [p, b])

    ### ### ### Adding a bond

    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)
    
    # Tests
    @assert p.length == origLength + 2  # Length increased by 2
    @assert haskey(p.pos, 17)  # New positions exist
    @assert haskey(p.pos, 18)
    @assert p.pos[5][5] == 17  # Check connections were updated
    @assert p.pos[17][5] == 18
    @assert p.pos[18][5] == 6
    @assert p.pos[18][4] == 17
    @assert p.pos[17][4] == 5
    @assert p.pos[5][4] == 4
    @assert p.latS[5,8,9] == 1  # New lattice sites are occupied
    @assert p.latS[6,8,9] == 1

    ### ### ### Kinking a bond

    p.moveIdx = 4
    identifyBond!(p)
    p.moveType = 2  # kinks
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)
    
    # Tests
    @assert p.length == origLength  # Length unchanged
    @assert p.pos[5][5] == 17  # Check connections remain the same
    @assert p.pos[17][5] == 18
    @assert p.pos[18][5] == 6
    @assert p.pos[18][4] == 17
    @assert p.pos[17][4] == 5
    @assert p.pos[5][4] == 4
    @assert p.latS[4,8,9] == 1  # New lattice site is occupied
    @assert p.latS[5,8,8] == 0  # Old lattice site is zero

    p.moveIdx = 17
    identifyBond!(p)
    p.moveType = 2  # kinks
    p.moveVec = [0,0,-1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)
    
    # Tests
    @assert p.length == origLength  # Length unchanged
    @assert p.pos[5][5] == 17  # Check connections remain the same
    @assert p.pos[17][5] == 18
    @assert p.pos[18][5] == 6
    @assert p.pos[18][4] == 17
    @assert p.pos[17][4] == 5
    @assert p.pos[5][4] == 4
    @assert p.latS[6,8,9] == 0  # New lattice site is occupied
    @assert p.latS[6,8,8] == 1  # Old lattice site is zero

    ### ### ### Removing a bond

    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 3  # remove
    p.moveVec = [0,0,-1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)
    
    # Tests
    @assert p.length == origLength - 2  # Length reduced by 2
    @assert p.pos[4][5] == 18  # Check connections updated
    @assert p.pos[18][5] == 6
    @assert p.pos[6][4] == 18
    @assert p.pos[18][4] == 4
    @assert p.latS[5,8,9] == 0  
    @assert p.latS[5,8,8] == 1

    # v = view(p.pos[18], 1:3)
    # p.latS[v...] = 0

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
    global p = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # Bulk monomers
    global b = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(p)
    initBulk!(b)
    initLattice!(l, [p, b])

    ### ### ### Adding a bond

    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)

    p.moveIdx = 3
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [1,1,1,1]

    p.moveIdx = 4
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [1,1,2,1]

    p.moveIdx = 5
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [2,2,1,1]

    p.moveIdx = 17
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [1,1,1,3]

    p.moveIdx = 18
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [2,2,1,1]

    
    p.moveIdx = 18
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,1,0]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    makeMove!(p)

    p.moveIdx = 18
    identifyBond!(p)
    moveTypes!(p)

    @assert p.moveTypes == [1,2,1,2]

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
    global p = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # Bulk monomers
    global b = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    # initializing lattice and species
    initPoly!(p)
    initBulk!(b)
    initLattice!(l, [p, b])

    ### ### ### Adding a bond

    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    makeMove!(p);

    p.moveIdx = 17
    identifyBond!(p)
    moveTypes!(p)

    @show p.moveTypes

    proposeMove!(p)

    @show p.moveType

    println("")

end

## 1010011111110000100000101010111100010100001100101110000
# Testing move
@time let
    date = "20251126"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    a = Bulk(;
        length = 27,
        intPot = -1,
        intPotP = -1,
        chemPot = 1,
        initShape = "box",
        initParams = 1,
        initPos = [8,8,8]
    )

    b = Bulk(;
        length = 32,
        intPot = -1,
        intPotP = -1,
        chemPot = 1,
        initShape = "uniform"
    )

    b.partner = a
    a.partner = b

    slist = [a, b]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0]
    frame = [1,1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0,0]
    niter = 0
    moveMade = [false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hams = [hamsA, hamsB]
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaH = [accDeltaHA, accDeltaHB]
    
    # first frame, initial energy
    push!(hams[1], a.ham)
    push!(hams[2], b.ham)
    push!(accDeltaH[1], a.ham)
    push!(accDeltaH[2], b.ham)

    ham!(a)
    ham!(b)
    a.deltaH = 0
    b.deltaH = 0

    a.name = string(a.name, "_a")
    b.name = string(b.name, "_b")
    color = ["blue (A)", "orange (B)"]

    # println("frame ", frame_joint, ", ", color[1], " init (", frame[1], ")")
    # println("frame ", frame_joint, ", ", color[2], " init (", frame[2], ")")
    exportSpecies!(slist, string(outpath, "_frame$(frame_joint)"))

    while frame_joint <= nframes - 2
        for (i, s) in enumerate(slist)
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH
                ip = i == 1 ? 2 : 1

                # update deltaH for the partner
                deltaH[ip] += s.partner.deltaHP
                
                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], s.ham)

                frame[i] += 1
                frame_joint += 1

                # println("frame ", frame_joint, ", ", color[i], " moved (", frame[i], ")")
                # exportSpecies!(slist, string(outpath, "_frame$(frame_joint)"))
            end
        end
    end

    # gyration radii of a straight line and a sphere of volume 'length'
    exportSpecies!(slist, string(outpath, "_frame$(2)"))

    println("A acceptance rate: ", nmoves[1] / nframes)
    println("B acceptance rate: ", nmoves[2] / nframes)

    @assert hams[1] == accDeltaH[1]
    @assert hams[2] == accDeltaH[2]

    p1 = plot(eachindex(hams[1]), hams[1][:];
        label = "ham A",
        xlabel = "frame",
        ylabel = "Energy",
        color = "blue",
        legend = :topright)
    
    p2 = plot(eachindex(hams[2]), hams[2][:];
        label = "ham B",
        xlabel = "frame",
        ylabel = "Energy",
        color = "orange",
        legend = :topright)

    plot!(p1, eachindex(accDeltaH[1]), accDeltaH[1][:];
        label = "accDeltaH A",
        color = "blue",
        linestyle = :dash)

    plot!(p2, eachindex(accDeltaH[2]), accDeltaH[2][:];
        label = "accDeltaH B",
        color = "orange",
        linestyle = :dash)

    plot(p1, p2; layout = (2, 1), size = (900, 600))
end

## 1010011111110000100000101010111100010100001100101110000
# Testing hamiltonian functions
let
    date = "20250924"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 1
    )

    # polymer
    global p = Species(;
        length = 16,
        maxlength = 32,
        boxlength = l.boxlength,
        intPot = 1.0,
        chemPot = 1.0
    )

    # initializing lattice and species
    initPoly!(p)
    initLattice!(l, [p])

    # initialize hamiltonians
    initHamPoly!(p)
    hamInitial = p.ham
    @show hamInitial
    @show p.chemPot, p.intPot

    ### ### ### Adding a bond
    
    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    deltaHamPoly!(p)
    @show p.deltaH
    @show hamInitial + p.deltaH
    
    makeMove!(p)
    p.length = length(keys(p.pos))

    hamPoly!(p)
    @show p.ham
    @assert p.ham == hamInitial + p.deltaH

    # ### ### ### Removing a bond

    hamInitial = p.ham

    p.moveIdx = 17
    identifyBond!(p)
    p.moveType = 3  # bond addition
    p.moveVec = [0,0,-1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec

    deltaHamPoly!(p)
    @show p.deltaH
    @show hamInitial + p.deltaH

    makeMove!(p)
    p.length = length(keys(p.pos))

    hamPoly!(p)
    @show p.ham
    @assert p.ham == hamInitial + p.deltaH

    export!(l, outpath)
end

## 1010011111110000100000101010111100010100001100101110000
# Testing center of mass

let 
    date = "20250916"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 32,
        nspecies = 2
    )

    # polymer
    global p = Species(;
        length = 32,
        maxlength = 300,
        boxlength = l.boxlength,
        intPot = -2,
        chemPot = 5
    )

    # Bulk monomers
    global b = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength,
        intPot = 0.0,
        chemPot = 0.0
    )

    # initializing lattice and species
    initPoly!(p)
    initBulk!(b)
    initLattice!(l, [p, b])

    @show centerMass(p) 

    export!(l, outpath)

end

## 1010011111110000100000101010111100010100001100101110000
# Multiple method definitions for a function

let 
    abstract type Person end
    abstract type PersonA <: Person end
    abstract type PersonB <: Person end

    mutable struct A <: PersonA
        name::String
        ssn::Int
    end

    mutable struct B <: PersonB
        name::String
        ssn::Int
    end

    mutable struct C <: Person
        name::String
    end

    function f(obj::A)
        println(obj.name, obj.ssn)
    end

    function f(obj::B)
        println(obj.name, obj.ssn)
    end

    function f(obj::C)
        println(obj.name)
    end

    function g(obj::Person)
        println(obj.ssn)
    end

    personA = A("Alice", 1)
    personB = B("Bob", 2)
    personC = C("Clarence")
    g(personB)

    @time hasproperty(personC, :ssn)
end

## 1010011111110000100000101010111100010100001100101110000
# Multi-component interactions

# to test:
# ~ bulk-bulk interaction: hamiltonian
# ~ fastbulk: movableIdxs
# * fastbulk: hastingsFactor

# to implement
# * bulk: static species
# * bulk: initialize surface
# * bulk-poly interaction
# * fastPoly
# * kink-only polymer moves


@time let
    date = "20251202"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)



    exportSpecies!(slist, string(outpath, "_frame$(2)"))
end