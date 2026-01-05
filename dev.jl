let 
    using Droplet
    using StatsBase
    using Distributions
    using Plots
    using Profile
    using OrderedCollections
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
        boxlength = 4,
        nspecies = 2
    )

    a = FastBulkMonomer(;
        length = 4,
        initShape = "box",
        initParams = 1,
        initPos = [2,2,2]
    )

    b = FastBulkMonomer(;
        length = 2,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 2,
        initShape = "uniform"
    )

    # self and interaction potentials
    intPot!(a, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)
    intPot!(a, b, -1.0)
    intPot!(a, c, -1.0)

    slist = [a, b, c]
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
    push!(hams_[1], a.ham)
    push!(hams_[2], b.ham)
    push!(hams_[3], c.ham)
    push!(accDeltaH_[1], a.ham)
    push!(accDeltaH_[2], b.ham)
    push!(accDeltaH_[3], c.ham)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])
    end

    ham!(a)
    ham!(b)
    ham!(c)
    a.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0

    color = ["blue (A)", "orange (B)", "green (C)"]
    exportSpecies!(slist, string(outpath, "_frame$(1)"))

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH

                # update deltaH for the partners
                for (p, j) in s.intSpecies
                    if j != s.intIdx
                        deltaH[j] += p.deltaHP
                        # println(p.name, " ", deltaH[j])
                    end
                end

                # println(s.name, " ", deltaH[i])
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
    # exportSpecies!(slist, string(outpath, "_frame$(2)"))

    println("A acceptance rate: ", nmoves[1] / nframes)
    println("B acceptance rate: ", nmoves[2] / nframes)

    @assert hams[1] == accDeltaH[1]
    @assert hams[2] == accDeltaH[2]
    @assert hams[3] == accDeltaH[3]

    plots = []
    colors = ["blue", "orange", "green"]

    for i in eachindex(hams)
        pl = plot(eachindex(hams[i]), hams[i][:];
        label = "ham",
        xlabel = "frame",
        ylabel = "Energy",
        color = colors[i],
        legend = :topright)

        plot!(pl, eachindex(accDeltaH[i]), accDeltaH[i][:];
        label = "accDeltaH",
        color = colors[i],
        linestyle = :dash)

        push!(plots, pl)
    end

    plot(plots[1], plots[2]; layout = (2, 1), size = (600, 600))
end

## 1010011111110000100000101010111100010100001100101110000
# Function dispatch

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

    @time hasproperty(persontf, :ssn)
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
# * multi-species poly
#   * new type: subspecies, without p.pos and move-related containers
#   * keep a p.nspecies-vector of subspecies and a p.length-vector of types 
#   * access the subspecies at bead i in deltaHam when you need intPotP 
#   * other species use poly subspecies as partners
# * site occupation time-series and mean -> BP estimate of concentration
# * distribution of densities
#   * average histogram of lattice densities over all sites, should bifurcate
# * exp: nucleation sites prevent Ostwald ripening (insurance)
#   * plot n_enhancers vs. dilute density compared with critical concentration 

# * fastpoly (no longer relevant)
#   * store bead indices at lattice sites
#   * kink-only mode: movableIdxs: a bead can be kinked if 
#       - 1. it isn't on a corner
#       - 2. if it is one bond away from a corner
#       - 3. it has no nn in the direction of the corner
#     after a kink:
#       - kinked bead no longer movable
#       - its 4 new nns are no longer movable
#       - its neighbor along the chain becomes movable
#       - its neighbors' old nns become movable if they were stuck bcs of #3
#   * hastings factor
# * reproduce Mason's results


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
