### physical parameters

l = Lattice(;
    boxlength = 16,
    nspecies = 1
)

# bulk species
b = FastBulk(;
    length = 2048
)

# list containing all species
slist = [b]

# minimum size of a droplet for detection
minsize = 50

### simulation parameters

t_max = 10    # runtime in minutes
step = 1000

vals = range(0.1, stop=1.6, step=0.1)   # parallelized parameters (int. pot.)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [1.0]   # temperatures
betas = [0.0]    # (physical) chemical and interaction potentials
nparams = length(alphas)


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, temp, beta = args

    b = slist[1]
    intPot!(b, - 1.0 / temp)

    # initialize lattice and species
    init!(l, [b])

    nmoves = 0
    iter = 0
    frame = 1

    ndroplets = []
    droplet_sizes = []
    droplet_density = []
    dilute_density = []
    hams = []
    accs = []

    deltaH = 0
    moveMade = false

    droplets = filter(x -> length(x) >= minsize, cc(b))
    ddens = dropletDensities(b, droplets)

    push!(ndroplets, length(droplets))
    push!(droplet_sizes, map(x -> length(x), droplets))
    push!(droplet_density, mean(ddens))
    push!(dilute_density, diluteDensity(b, droplets))
    push!(hams, b.ham)   # initial energy
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10

        if moveMade
            deltaH += b.deltaH
        
            if mod(nmoves, step) == 0
                frame += 1

                ccs = cc(b)
                droplets = filter(x -> length(x) >= minsize, ccs)
                ddens = dropletDensities(b, droplets)
                
                push!(ndroplets, length(droplets))
                push!(droplet_sizes, map(x -> length(x), droplets))
                push!(droplet_density, mean(ddens))
                push!(dilute_density, diluteDensity(b, droplets))
                push!(hams, hams[end] + deltaH)
                push!(accs, iter == 0 ? 0 : nmoves / iter)
                
                deltaH = 0
            end
        end
        
        moveMade = move!(b)
        nmoves += moveMade
        iter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, b.length, 4)

    res["species_name"] = b.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = temp
    res["beta"] = beta
    res["boxlength"] = b.boxlength
    res["length"] = b.length
    res["ndroplets"] = ndroplets
    res["droplet_sizes"] = droplet_sizes
    res["droplet_density"] = droplet_density
    res["dilute_density"] = dilute_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(b, coords)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
