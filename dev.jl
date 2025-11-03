let 
    using Droplet
    using StatsBase
    using Distributions
    using Plots
    using Profile
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
# Testing poly move
@time let
    date = "20251028"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 42,
        nspecies = 2
    )

    # polymer
    global p = Poly(;
        length = 32,
        maxlength = 256,
        boxlength = l.boxlength,
        intPot = 1,
        chemPot = 1,
        initPos = [5,21,21]
    )

    # Bulk monomers
    global b = Bulk(;
        length = 16,
        boxlength = l.boxlength,
        intPot = 0.0,
        chemPot = 0.0
    )
    slist = [p]

    # initialize lattice and species
    init!(p)
    init!(b)
    initLattice!(l, slist)
    
    # initialize hamiltonians
    initHam!(p)
    hamInitial = p.ham

    #### iterating multiple moves
    nmoves = 0
    niter = 0
    frame = 1
    nframes = 100
    step = 50
    deltaH = 0
    
    moveMade = false
    radii = zeros(nframes)
    lengths = zeros(nframes)
    hams = zeros(nframes)
    accDeltaH = zeros(nframes)
    moves = ["add", "kink", "remove"]

    # println("Frame ", frame)
    # @show p.moveIdx, p.moveType, p.moveVec, p.length
    # @show p.pos
    
    # first frame
    radii[frame] = radius(p)
    lengths[frame] = p.length
    accDeltaH[frame] = hamInitial
    
    ham!(p)
    hams[frame] = hamInitial
    
    # @show p.ham, p.deltaH
    # println("")

    exportSpecies!(slist, string(outpath, "_frame$(frame)"))

    while frame <= nframes - 1
        
        if moveMade

            deltaH += p.deltaH

            if mod(nmoves, step) == 0
            
                frame += 1

                # println("Frame ", frame)
                # @show p.moveIdx, moves[p.moveType], p.moveVec, p.length
                # @show p.pos
                radii[frame] = radius(p)
                lengths[frame] = p.length
                accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                deltaH = 0
                
                ham!(p)
                hams[frame] = p.ham

                # @show p.ham, p.deltaH
                # println("")
                
                # export!(l, string(outpath, "_frame$(frame)"))
                exportSpecies!(slist, string(outpath, "_frame$(frame)"))

                if hams[frame] != accDeltaH[frame]
                    println("Frame ", frame)
                    @show p.moveIdx, moves[p.moveType], p.moveVec, p.length
                    @show p.pos
                    @show p.ham, p.deltaH
                    @assert hams[frame] == accDeltaH[frame]
                end
            end
        end
        
        moveMade = move!(p)
        nmoves += moveMade
        niter += 1
    end

    # gyration radii of a straight line and a sphere of volume 'length'
    rg_line = radiusLine(p.boxlength)
    rg_sphere = radiusSphere(length(keys(p.pos)))

    println("Acceptance rate: ", nmoves / niter)

    ham!(p)
    println("Hamiltonian: ", p.ham)

    # plot(1:nframes, radii, label = "radius of gyration")
    # plot!(1:nframes, lengths, label = "polymer length")

    p1 = plot(1:nframes, radii;
        label = "radius of gyration",
        xlabel = "frame",
        ylabel = "Rg",
        legend = :topright)

    p2 = plot(1:nframes, hams;
        label = "ham",
        xlabel = "frame",
        ylabel = "Energy",
        legend = :topright)
    plot!(p2, 1:nframes, accDeltaH;
        label = "accDeltaH",
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
# What is context?

let 
    filename = "ctx.json"
    open(filename, "w") do f
        JSON.print(f, CONTEXT, 4)
    end
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
