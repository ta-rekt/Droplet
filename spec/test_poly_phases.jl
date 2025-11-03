### physical parameters

global l = Lattice(;
    boxlength = 42,
    nspecies = 1
)

# polymer
global p = Poly(;
    length = 32,
    maxlength = 512,
    boxlength = boxlength,
    intPot = nothing,
    chemPot = nothing,
    initPos = [5,21,21]
)

# list containing all species
slist = [p]


### simulation parameters

nframes = 10000
step = 100

vals = range(-1, stop=1, step=0.25)   # parallelized parameters (int. pot.)
# alphas = zeros(9)
nrepeats = 1

# for (i, v) in enumerate(vals)
#     alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
# end

alphas = [-1.0]
betas = [1.0]    # chemical potential
nparams = length(alphas)


### functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    p = slist[1]
    p.intPot = alpha
    p.chemPot = beta

    # initialize lattice and species
    init!(p)
    initLattice!(l, [p])

    # initialize hamiltonians
    initHam!(p)

    nmoves = 0
    niter = 0
    frame = 1

    radii = zeros(nframes)
    lengths = zeros(nframes)
    hams = zeros(nframes)
    accs = zeros(nframes)
    deltaH = 0
    moveMade = false
    
    radii[frame] = radius(p)
    lengths[frame] = p.length
    hams[frame] = p.ham   # initial energy
    accs[frame] = 0
    
    # run simluation
    while frame <= nframes - 1

        if moveMade
            deltaH += p.deltaH
        
            if mod(nmoves, step) == 0
                
                frame += 1
                
                radii[frame] = radius(p)
                lengths[frame] = p.length
                hams[frame] = hams[frame-1] + deltaH
                accs[frame] = niter == 0 ? 0 : nmoves / niter
                
                deltaH = 0
            end
        end
        
        moveMade = move!(p)
        nmoves += moveMade
        niter += 1
    end

    # collect results
    coords = Matrix{Int32}(undef, p.length, 3)

    res["species_name"] = p.name
    res["nframes"] = nframes
    res["step"] = step
    res["alpha"] = alpha
    res["beta"] = beta
    res["maxlength"] = p.maxlength
    res["boxlength"] = p.boxlength
    res["gyration_radius"] = radii
    res["polymer_length"] = lengths
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(p, coords)
    res["runtime_s"] = (time_ns() - t_start) / 1e9

    return res
end
