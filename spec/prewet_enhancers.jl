### physical parameters

l = Lattice(;
    boxlength = 16,
    nspecies = 3
)

## species initial configuration

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

slist = [a1, a2, b, c]

### simulation parameters

nframes = 10000
step = 1000

# minimum size of a droplet for detection
minsize = 5

vals = range(124, stop=1024, step=60)  # parallelized parameters (int. pot.)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [128]   # number of particles
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)
nspecies = 4
nmoving = 1

### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    a1, a2, b, c = slist

    c.length = alpha

    # self-interactions
    intPot!(a1, -1.0)
    intPot!(a2, -1.0)
    intPot!(b, -1.0)
    intPot!(c, -1.0)

    # multi-species interactions
    intPot!(c, a1, -0.5) 
    intPot!(c, a2, -0.5)   
    intPot!(c, b, -2.0) 

    slist = [a1, a2, b, c]         # all species
    movingsps = [(c, c.intIdx)]    # species that move

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = 0
    frame = ones(Int, nmoving)
    frame_joint = 1
    deltaH = zeros(Float64, nspecies)

    ndroplets = []
    droplet_sizes = []
    droplet_densities = []
    dilute_density = []
    hams = []
    accs = []

    niter = 0
    moveMade = falses(nmoving)

    for (sp, idx) in a1.intSpecies
        hams_sp = Vector{Float64}(undef, 0)

        push!(hams_sp, sp.ham)
        push!(hams, hams_sp)

        ham!(sp)
        sp.deltaH = 0.0
    end

    while frame_joint <= nframes - 2
        for (s, i) in movingsps
            if  move!(s)
                nmoves += 1

                if mod(nmoves, step) == 0
                    deltaH[i] += s.deltaH

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end
                    
                    ccs = cc(b)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    
                    push!(ndroplets, length(droplets))
                    push!(droplet_sizes, map(x -> length(x), droplets))
                    push!(droplet_densities, dropletDensities(b, droplets))
                    push!(dilute_density, diluteDensity(b, droplets))
                    push!(accs, niter == 0 ? 0 : nmoves / niter)

                    push!(hams[i], hams[i][end] + deltaH[i])
                    deltaH[i] = 0

                    frame[i] += 1
                    frame_joint += 1
                end
            end
        end
    end

    # collect results
    coords = Matrix{Int32}(undef, c.length, 4)

    res["species_name"] = c.name
    res["nframes"] = nframes
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["boxlength"] = c.boxlength
    res["length"] = c.length
    res["ndroplets"] = ndroplets
    res["droplet_sizes"] = droplet_sizes
    res["droplet_densities"] = droplet_densities
    res["dilute_density"] = dilute_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(c, coords)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end