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
        boxlength = 16,
        nspecies = 2
    )

    c = FastBulk(;
        length = 256,
        initShape = "uniform"
    )

    a1 = FastBulk(;
        length = 5,
        initShape = "line",
        initParams = 1,
        initPos = [1,5,8]
    )

    a2 = FastBulk(;
        length = 5,
        initShape = "line",
        initParams = 1,
        initPos = [12,5,8]
    )

    b = FastBulk(;
        length = 6,
        initShape = "line",
        initParams = 1,
        initPos = [6,5,8]
    )

    # self-interaction
    intPot!(a1, -1.0)
    intPot!(a2, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)

    # inter-species interaction
    intPot!(c, a1, -0.5)  # a-c interaction
    intPot!(c, a2, -0.5)  # b-c interaction  
    intPot!(c, b, -2.0)  # a-b interaction (optional)

    slist = [a1, a2, b, c]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = zeros(Int32, 4)
    frame = ones(Int, 4)
    frame_joint = 1
    nframes = 1000
    step = 1
    deltaH = zeros(Float64, 4)
    niter = 0
    moveMade = falses(4)

    hams = []
    accDeltaH = []

    for (sp, idx) in a1.intSpecies
        hams_sp = Vector{Float64}(undef, 0)
        accDeltaH_sp = Vector{Float64}(undef, 0)

        push!(hams_sp, sp.ham)
        push!(accDeltaH_sp, sp.ham)

        push!(hams, hams_sp)
        push!(accDeltaH, accDeltaH_sp)

        ham!(sp)
        sp.deltaH = 0.0
    end

    while frame_joint <= nframes - 2
        for (s, i) in a1.intSpecies
            if s.intIdx == c.intIdx && move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH

                # update deltaH for the partners
                for (p, j) in s.intSpecies
                    if j != s.intIdx
                        deltaH[j] += p.deltaHP
                    end
                end

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], s.ham)

                frame[i] += 1
                frame_joint += 1
            end
        end
    end

    # gyration radii of a straight line and a sphere of volume 'length'
    exportSpecies!(slist, string(outpath, "_frame$(2)"))

    println("A acceptance rate: ", nmoves[1] / nframes)
    println("B acceptance rate: ", nmoves[2] / nframes)

    @assert hams[1] == accDeltaH[1]
    @assert hams[2] == accDeltaH[2]

    plots = []
    colors = ["blue", "red", "orange", "green"]

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

    plot(plots[1]; layout = (1, 1), size = (300, 600))
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

    @time hasproperty(personC, :ssn)
end

## 1010011111110000100000101010111100010100001100101110000
# Multi-component interactions

# to test:
# ~ intPot! function
#   ~ give init a list of interacting species, int. matrix (symmetric)
#   ~ each species is given a copy of the matrix, partner species and idxs
#   ~ partner intPots are taken from the matrix
#   ~ update hamiltonian functions
# ~ multi-species interactions
# * fastbulk flory huggins (insurance)

# to implement

# * filament prewetting (insurance)
# * effective force between two prewet enhancers (insurance)
# * nucleation sites prevent Ostwald ripening (insurance)
# monday 08/12
# * poly: length-preserving moves
#   * adapt the proposal probabilities
#   * run a burn-in with ADD and REMOVE, then equilibrate 
#     with KINK only once length = lmax. rebuild the pos dict with
#     keys in 1:lmax
# * fastpoly
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
# * bulk-poly interaction (uniform)
# * multi-species poly
#   * new type: subspecies, without p.pos and move-related containers
#   * keep a p.nspecies-vector of subspecies and a p.length-vector of types 
#   * access the subspecies at bead i in deltaHam when you need intPotP 
#   * other species use poly subspecies as partners
# monday 15/12
# * reproduce Mason's results
# thursday 18/12


@time let
    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    ############ droplet boundary ############

    l = Lattice(;
        boxlength = 16,
        nspecies = 3
    )

    c = FastBulk(;
        length = 128,
        initShape = "uniform"
    )

    a1 = FastBulk(;
        length = 5,
        initShape = "line",
        initParams = 1,
        initPos = [1,5,8]
    )

    a2 = FastBulk(;
        length = 5,
        initShape = "line",
        initParams = 1,
        initPos = [12,5,8]
    )

    b = FastBulk(;
        length = 6,
        initShape = "line",
        initParams = 1,
        initPos = [6,5,8]
    )

    intPot!(a1, -1.0)
    intPot!(a2, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)

    # Set up multi-species interactions using the new intPot! function
    intPot!(c, a1, -0.5)  # a-c interaction
    intPot!(c, a2, -0.5)  # b-c interaction  
    intPot!(c, b, -2.0)  # a-b interaction (optional)

    slist = [a1, a2, b, c]
    
    # initialize lattice and species
    init!(l, slist)

    println(map(x -> x[2], a1.intSpecies))

    exportSpecies!(slist, string(outpath, "_frame$(1)"))

end