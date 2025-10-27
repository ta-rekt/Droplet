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
    outpath = "data/$(outdir)/"
    Base.Filesystem.mkpath(outpath)

    ### specification file
    specpath = string("../spec/", specfile)
    include(specpath)

    nparams_per_cpu = div(nparams, ncpus)

    ### output data
    results = Dict{Int64, Dict{Int64, Dict{String, Vector{Float32}}}}()

    for k in 1:nparams_per_cpu

        keya = (cpukey - 1) * nparams_per_cpu + k
        alpha = alphas[keya]
        results[keya] = Dict{Int64, Dict{String, Vector{Float32}}}()

        for (keyb, beta) in enumerate(betas)

            res = Dict{String, Vector{Float32}}()
            args = [slist, alpha, beta]

            results[keya][keyb] = target!(res, args...)
        end
    end

    
    ### write data
    file = string(outpath, "sweep_output$(cpukey).json")

    open(file, "w") do f
        JSON.print(f, results, 4)
    end
end
