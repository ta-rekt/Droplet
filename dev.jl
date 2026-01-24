let 
    using Droplet
    using StatsBase
    using Distributions
    using Plots
    using Profile
    using OrderedCollections
end


## 1010011111110000100000101010111100010100001100101110000
# Testing Hamiltonian - PolyMonomer implementation
@time let
    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 4,
        nspecies = 1
    )

    p = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 1,
        initPos = [1,2,2],
        pattern=:uniform
    )
    
    initMonomers!(p)
    a = p.monomers[1]

    # self and interaction potentials
    intPot!(a, -1.0)

    slist = [p]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0]
    frame = [1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0.0]
    niter = 0
    moveMade = [false]

    hamsA = Vector{Float64}(undef, 0)
    hams_ = [hamsA]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], getHam(a))
    push!(accDeltaH_[1], getHam(a))

    bulkSpecies = Vector{Any}(undef, 0)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])

        if typeof(sp) != PolyMonomer
            push!(bulkSpecies, (sp, idx))
        end
    end

    ham!(a)
    p.deltaH = 0.0

    color = ["blue (A)"]
    exportSpecies!(slist, string(outpath, "_frame$(1)"))

    movetypes = ["add", "kink", "remove"]
    bad_frame = 1

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies

            s_prev = deepcopy(p)
            slist_prev = [s_prev]

            if move!(s)
                nmoves[i] += 1
                deltaH[i] += getDeltaHam(s)

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], getHam(s))

                frame[i] += 1
                frame_joint += 1

                if hams[1] != accDeltaH[1]
                    bad_frame = frame_joint
                    println("Inconsistency at frame: ", frame_joint)
                    println("Move type: ", movetypes[p.moveType])
                    println("Bead: ", p.moveIdx)

                    exportSpecies!(slist_prev, string(outpath, "_frame$(2)"))
                    exportSpecies!(slist, string(outpath, "_frame$(3)"))
                    frame_joint = nframes
                end
            end
        end
    end

    # plots = []
    # colors = ["blue"]
    # hamy = hams[1][bad_frame-1:bad_frame+9]
    # accy = accDeltaH[1][bad_frame-1:bad_frame+9]

    # pl = plot(eachindex(hamy), hamy;
    # label = "ham",
    # xlabel = "frame",
    # ylabel = "Energy",
    # color = colors[1],
    # legend = :topright)

    # plot!(pl, eachindex(accy), accy;
    # label = "accDeltaH",
    # color = colors[1],
    # linestyle = :dash)

    # push!(plots, pl)

    # plot(plots[1]; layout = (1, 1), size = (600, 300))
end

## 1010011111110000100000101010111100010100001100101110000
# Testing Hamiltonian and its increment
@time let
    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 16,
        nspecies = 3
    )

    p = Polymer(;
        length = 16,
        maxlength = 64,
        nMonomerTypes = 1,
        initPos = [1,8,8],
        pattern=:uniform
    )
    
    initMonomers!(p)
    a = p.monomers[1]

    b = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    # self and interaction potentials
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -1.0)
    intPot!(a, c, -1.0)

    slist = [p, b, c]
    # slist = [a, b]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0,0]
    frame = [1,1,1]
    frame_joint = 1
    nframes = 1000
    step = 1
    deltaH = [0.0,0.0,0.0]
    niter = 0
    moveMade = [false, false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hamsC = Vector{Float64}(undef, 0)
    hams_ = [hamsA, hamsB, hamsC]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaHC = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], getHam(a))
    push!(hams_[2], getHam(b))
    push!(hams_[3], getHam(c))
    push!(accDeltaH_[1], getHam(a))
    push!(accDeltaH_[2], getHam(b))
    push!(accDeltaH_[3], getHam(c))

    bulkSpecies = Vector{Any}(undef, 0)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])

        if typeof(sp) != PolyMonomer
            push!(bulkSpecies, (sp, idx))
        end
    end

    ham!(a)
    ham!(b)
    ham!(c)
    p.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0
    p.deltaHP = 0.0
    b.deltaHP = 0.0
    c.deltaHP = 0.0

    color = ["blue (A)", "orange (B)", "green (C)"]
    exportSpecies!(slist, string(outpath, "_frame$(1)"))

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += getDeltaHam(s)

                # update deltaH for the partners
                for (partner, j) in a.intSpecies
                    if j != i
                        deltaH[j] += getDeltaHP(partner)
                        # @show p.name, deltaH[j]
                    end
                end

                # @show s.name, deltaH[i] 
                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], getHam(s))

                frame[i] += 1
                frame_joint += 1

                # println("frame ", frame_joint, ", ", color[i], " moved (", frame[i], ")")
                # exportSpecies!(slist, string(outpath, "_frame$(frame_joint)"))

                if hams[1] != accDeltaH[1] || hams[2] != accDeltaH[2] || hams[3] != accDeltaH[3]
                    println("inconsistency at frame ", frame_joint)
                    break
                end
            end
        end
    end

    # gyration radii of a straight line and a sphere of volume 'length'
    # exportSpecies!(slist, string(outpath, "_frame$(2)"))

    println("A acceptance rate: ", nmoves[1] / nframes)
    println("B acceptance rate: ", nmoves[2] / nframes)

    @assert hams[1] == accDeltaH[1]
    @assert hams[2] == accDeltaH[2]
    @assert hams[3] == accDeltaH[3]

    # plots = []
    # colors = ["blue", "orange", "green"]

    # for i in eachindex(hams)
    #     pl = plot(eachindex(hams[i]), hams[i][:];
    #     label = "ham",
    #     xlabel = "frame",
    #     ylabel = "Energy",
    #     color = colors[i],
    #     legend = :topright)

    #     plot!(pl, eachindex(accDeltaH[i]), accDeltaH[i][:];
    #     label = "accDeltaH",
    #     color = colors[i],
    #     linestyle = :dash)

    #     push!(plots, pl)
    # end

    # plot(plots[1], plots[2], plots[3]; layout = (3, 1), size = (600, 600))
end

## 1010011111110000100000101010111100010100001100101110000
# Multi-component interactions

# to test:
# * exp: fastbulkMonomer flory huggins (insurance)
# * exp: filament prewetting (insurance)
#   * vary concentration and look at axial/radial density, compare with control
# * exp: prewet enhancers (insurance)
#   * radial densities of two nearby enhancers
#   * effective force, comparison with electrostatic

# to implement
# * site occupation time-series and mean -> BP estimate of concentration
# * distribution of densities
#   * average histogram of lattice densities over all sites, should bifurcate
# * exp: nucleation sites prevent Ostwald ripening (insurance)
#   * plot n_enhancers vs. dilute density compared with critical concentration 

# * reproduce Mason's results

# TODO:
# the idea: I'd like to be able to parallelize setting up experiments and 
# working on new features. I need to develop a small preliminary workflow
# where I can describe in detail the experiment I need to set up, and have 
# the agent make a proposal. the experiment proposal format needs to be
# designed, it must consist at mimimum of: creating a spec file, a summary
# of the experiment, and a description of the resulting plots, expected
# results, the technical specs such as number of cores to run on, expected
# simulation time, varied parameters, type of output data, etc.

# - set up paper airplane agentic experimental pipeline (run locally for now)
#   - run design iterations to produce description files of the experiment 
#     design/execution pipeline. write detailed context, let it ask questions
# - try it on a dummy experiment for validation
# - use it to design the bulk-mediated polymer collapse experiment
# - run bulk-mediated polymer collapse experiment
# - implement heterogeneous polymer
# - test heterogeneous polymer
# - expand agentic pipeline to cluster


@time let
    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    ############ droplet boundary ############

    l = Lattice(;
        boxlength = 20,
        nspecies = 3
    )

    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )

    # rnap = FastBulkMonomer(;
    #     length = 32,
    #     initShape = "uniform"
    # )

    ### rod 1
    dna1 = FastBulkMonomer(;
        length = 7,
        initShape = "line",
        initParams = 1,
        initPos = [1,5,8]
    )

    enhancer = FastBulkMonomer(;
        length = 6,
        initShape = "line",
        initParams = 1,
        initPos = [8,5,8]
    )

    dna2 = FastBulkMonomer(;
        length = 7,
        initShape = "line",
        initParams = 1,
        initPos = [14,5,8]
    )


    ### rod 1
    dna1_l = FastBulkMonomer(;
        length = 7,
        initShape = "line",
        initParams = 1,
        initPos = [1,10,8]
    )

    enhancer_l = FastBulkMonomer(;
        length = 6,
        initShape = "line",
        initParams = 1,
        initPos = [8,10,8]
    )

    dna2_l = FastBulkMonomer(;
        length = 7,
        initShape = "line",
        initParams = 1,
        initPos = [14,10,8]
    )

    # # self-interactions
    # intPot!(tf, -1.0)
    # intPot!(rnap, -1.0)

    # # multi-species interactions
    # intPot!(tf, dna1, -0.5)  
    # intPot!(tf, dna2, -0.5)  
    # intPot!(tf, enhancer, -2.0)

    # intPot!(rnap, dna1, -0.5)  
    # intPot!(rnap, dna2, -0.5)  
    # intPot!(rnap, enhancer, -0.5)

    # intPot!(tf, rnap, -2.0)

    slist = [dna1, dna2, enhancer, dna1_l, dna2_l, enhancer_l, tf]
    
    # initialize lattice and species
    init!(l, slist)

    exportSpecies!(slist, string(outpath, "_frame$(1)"))
end

## 1010011111110000100000101010111100010100001100101110000
# Multi-component polymer

@time let

    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    ############ droplet boundary ############

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

    intPot!(tf, dna, 0.3)
    intPot!(tf, promoter, 2.0)

    intPot!(tf, 1.0)
    intPot!(dna, 0.5)
    intPot!(promoter, 0.5)

    # map(p -> (p.name, p.intPot, getIntPot(tf, p)), poly.monomers)

    # add all species here, including polymonomers
    # let init initialize them
    # let ham use a mono's parent poly's position, or a bulk mono's position
    species = [tf, dna, promoter]

    init!(l, species)

    @show getIntPot(tf, dna)
    @show getIntPot(tf, promoter)
    @show getIntPot(tf, tf)

end
