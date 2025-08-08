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
# Adding species to a lattice, moving individual bonds
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
    global p1 = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    # dilute monomers
    p2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    initPoly!(p1)
    initDilute!(p2)
    initLattice!(l, [p1, p2])
    
    p1.moveIdx = 15
    p1.moveType = 2  # bond addition
    p1.moveVec = [0,0,1]
    makeMove!(l, p1)

    p1.moveIdx = 17
    p1.moveType = 2
    p1.moveVec = [0,0,1]
    makeMove!(l, p1)

    p1.moveIdx = 15
    p1.moveType = 3
    p1.moveVec = [-1,0,0]
    _, _, me = makeMove!(l, p1)
    println("move executed? ", me)
    
    ## turn off a site
    # v = view(p1.pos[18], 1:3)
    # latS = view(l.lat, l.speciesDir[p1.id], :, :, :)
    # latS[v...] = 0

    export!(l, outpath)

end


## 1010011111110000100000101010111100010100001100101110000
# Testing moveTypes! and proposeMove!
let 
    date = "20250807"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    global l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # polymer
    global s = Species(;
        length = 16,
        maxlength = 20,
        boxlength = l.boxlength
    )

    initPoly!(s)
    initLattice!(l, [s])

    s.moveIdx = 4
    s.moveType = 1
    s.moveVec = [0,0,-1]
    makeMove!(l, s)

    s.moveIdx = 17
    s.moveVec = [0,0,0]



    # turn off a site
    # v = view(s.pos[17], 1:3)
    # latS = view(l.lat, l.speciesDir[s.id], :, :, :)
    # latS[v...] = 0

    @show findall(iszero, s.bondVec)

    moveTypes!(s, bi, bj, bprev, bnext)

    @show s.moveIdx
    @show bi, bj, bprev, bnext
    @show s.moveTypes

    export!(l, outpath)

end