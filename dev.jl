using Droplet

## 1010011111110000100000101010111100010100001100101110000
# Saving and loading
let 

    date = "20250718"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 64,
        nspecies = 3
    )

    initLatticeRandom!(l)
    export!(l, outpath)
    save!(l, outpath)

    l0 = Lattice(;
        boxlength = 64,
        nspecies = 3
    )

    load!(l0, outpath)
    # @assert l.lat == l0.lat

end