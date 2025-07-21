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
# Adding species to a lattice
let 
    date = "20250721"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "testpoly"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    p1 = Poly(;
        length = 16,
        maxlength = 64,
        boxlength = l.boxlength,
        init_pos = [0, 5, 14]
    )

    p2 = Poly(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength
    )

    initPoly!(p1)
    initDilute!(p2)
    initLattice!(l, [p1, p2])

    export!(l, outpath)

end