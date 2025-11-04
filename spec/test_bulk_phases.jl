### physical parameters

global l = Lattice(;
    boxlength = 16,
    nspecies = 1
)

# bulk species
global b = Bulk(;
    length = 512,
    intPot = nothing,
    chemPot = nothing
)

# list containing all species
slist = [b]


### simulation parameters

nframes = 10000  
step = 10000     

vals = range(-1, stop=1, step=0.25)   # parallelized parameters (int. pot.)
alphas = zeros(9)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [-1.0]
betas = [1.0]    # chemical potential
nparams = length(alphas)


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    b = slist[1]
    b.intPot = alpha
    b.chemPot = beta

    # initialize lattice and species
    init!(l, [b])

    nmoves = 0
    niter = 0
    frame = 1

    radii = zeros(nframes)
    hams = zeros(nframes)
    accs = zeros(nframes)
    deltaH = 0
    moveMade = false
    
    radii[frame] = radius(b)
    hams[frame] = b.ham   # initial energy
    accs[frame] = 0
    
    # run simluation
    while frame <= nframes - 1

        if moveMade
            deltaH += b.deltaH
        
            if mod(nmoves, step) == 0
                frame += 1
                
                radii[frame] = radius(b)
                hams[frame] = hams[frame-1] + deltaH
                accs[frame] = niter == 0 ? 0 : nmoves / niter
                
                deltaH = 0
            end
        end
        
        moveMade = move!(b)
        nmoves += moveMade
        niter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, b.length, 3)

    res["species_name"] = b.name
    res["nframes"] = nframes
    res["step"] = step
    res["alpha"] = alpha
    res["beta"] = beta
    res["boxlength"] = b.boxlength
    res["length"] = b.length
    res["gyration_radius"] = radii
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(b, coords)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end
