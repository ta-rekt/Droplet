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
# Adding multiple species to a lattice
let 
    date = "20250724"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    l = Lattice(;
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
    
    p1.moveIdx = 4
    p1.moveType = 2  # bond addition
    p1.moveVec = [0,0,1]
    
    initLattice!(l, [p1, p2])
    makeMove!(l, p1)

    p1.moveIdx = 17
    p1.moveType = 1  # removal
    p1.moveVec = [0,0,-1]
    makeMove!(l, p1)
    
    ## turn off a site
    # v = view(p1.pos[1], 1:3)
    # latS = view(l.lat, l.speciesDir[p1.id], :, :, :)
    # latS[v...] = 0

    export!(l, outpath)

end
