import Pkg
using JSON
using Statistics
using Droplet

let
    ### command-line arguments
    specfile = ARGS[1]
    cpukey = parse(Int32, ARGS[2])
    ncpus = parse(Int32, ARGS[3])
    outdir = ARGS[4]

    ### output directory
    Base.Filesystem.mkpath(outdir)

    ### specification file
    specpath = string("../", specfile)
    include(specpath)

    nparams_per_cpu = div(nparams, ncpus)

    ### output data
    results = Dict{Int64, Dict{Int64, Dict{String, Any}}}()
    file = string(outdir, "_cpu$(cpukey).json")

    for k in 1:nparams_per_cpu

        keya = (cpukey - 1) * nparams_per_cpu + k
        alpha = alphas[keya]
        results[keya] = Dict{Int64, Dict{String, Any}}()

        for (keyb, beta) in enumerate(betas)

            results[keya][keyb] = Dict{String, Any}()
            args = [slist, keya, alpha, keyb, beta]

            results = target!(results, args...)
        end
    end
end
