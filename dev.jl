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
    date = "20251202"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "test"
    outpath = string(outdir, outfile)

    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    a = FastBulk(;
        length = 16,
        intPot = -1,
        intPotP = -2,
        chemPot = 1,
        initShape = "line",
        initParams = 1,
        initPos = [1,8,8]
    )

    b = FastBulk(;
        length = 128,
        intPot = -1,
        intPotP = a.intPotP,
        chemPot = 1,
        initShape = "uniform"
    )

    b.partner = a
    a.partner = b

    slist = [a, b]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0]
    frame = [1,1]
    frame_joint = 1
    nframes = 1000
    step = 1
    deltaH = [0,0]
    niter = 0
    moveMade = [false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hams = [hamsA, hamsB]
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaH = [accDeltaHA, accDeltaHB]
    
    # first frame, initial energy
    push!(hams[1], a.ham)
    push!(hams[2], b.ham)
    push!(accDeltaH[1], a.ham)
    push!(accDeltaH[2], b.ham)

    ham!(a)
    ham!(b)
    a.deltaH = 0
    b.deltaH = 0

    a.name = string(a.name, "_a")
    b.name = string(b.name, "_b")
    color = ["blue (A)", "orange (B)"]

    # println("frame ", frame_joint, ", ", color[1], " init (", frame[1], ")")
    # println("frame ", frame_joint, ", ", color[2], " init (", frame[2], ")")
    exportSpecies!(slist, string(outpath, "_frame$(frame_joint)"))

    while frame_joint <= nframes - 2
        s = b
        i = 2
        if move!(s)
            nmoves[i] += 1
            deltaH[i] += s.deltaH
            ip = i == 1 ? 2 : 1

            # update deltaH for the partner
            deltaH[ip] += s.partner.deltaHP
            
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

    # gyration radii of a straight line and a sphere of volume 'length'
    exportSpecies!(slist, string(outpath, "_frame$(2)"))

    println("A acceptance rate: ", nmoves[1] / nframes)
    println("B acceptance rate: ", nmoves[2] / nframes)

    @assert hams[1] == accDeltaH[1]
    @assert hams[2] == accDeltaH[2]

    p1 = plot(eachindex(hams[1]), hams[1][:];
        label = "ham A",
        xlabel = "frame",
        ylabel = "Energy",
        color = "blue",
        legend = :topright)
    
    p2 = plot(eachindex(hams[2]), hams[2][:];
        label = "ham B",
        xlabel = "frame",
        ylabel = "Energy",
        color = "orange",
        legend = :topright)

    plot!(p1, eachindex(accDeltaH[1]), accDeltaH[1][:];
        label = "accDeltaH A",
        color = "blue",
        linestyle = :dash)

    plot!(p2, eachindex(accDeltaH[2]), accDeltaH[2][:];
        label = "accDeltaH B",
        color = "orange",
        linestyle = :dash)

    plot(p1, p2; layout = (2, 1), size = (900, 600))
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
# * intPot! function

# to implement

# * k-way species interactions
#   ~ give init a list of interacting species, int. matrix (symmetric)
#   ~ each species is given a copy of the matrix, partner species and idxs
#   ~ partner intPots are taken from the matrix
#   * update hamiltonian functions
# * prewetting experiments with rigid rods (insurance)
# * multicomponent phase separation, try mean-field (insurance)
# * nucleation sites prevent Ostwald ripening (insurance)
# * fastbulk speed test (insurance)
# * fastbulk droplet density (insurance)
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
    date = "20251203"
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
        length = 7,
        initShape = "uniform"
    )

    a = FastBulk(;
        length = 16,
        initShape = "line",
        initParams = 1,
        initPos = [1,5,8]
    )

    b = FastBulk(;
        length = 16,
        initShape = "line",
        initParams = 1,
        initPos = [1,12,8]
    )

    intPot!(a, -1.0)

    # Set up multi-species interactions using the new intPot! function
    intPot!(a, c, -3.0)  # a-c interaction
    intPot!(b, c, -2.0)  # b-c interaction  
    intPot!(a, b, -0.5)  # a-b interaction (optional)

    slist = [a, b, c]
    
    # initialize lattice and species
    init!(l, slist)

    println(map(x -> x[2], a.intSpecies))

    exportSpecies!(slist, string(outpath, "_frame$(1)"))

end