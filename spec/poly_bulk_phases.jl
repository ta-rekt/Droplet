### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 16,
    nspecies = 2
)

# polymer species
poly = Polymer(;
    length = 16,
    maxlength = 64,
    initShape = "line",
    initPos = [8, 8, 8]
)

# bulk species
bulk = FastBulk(;
    length = 1024,
    initShape = "uniform"
)

# list containing all species
slist = [poly, bulk]


### simulation parameters
# MARK: -Sim. params

nframes = 50000
step = 10000     

vals = range(-1, stop=1, step=0.25)   # parallelized parameters (polymer self-interaction)
alphas = zeros(9)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

# alphas = [-1.0]
betas = [0.0]    # placeholder
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, alpha, beta = args

    poly, bulk = slist
    
    # set interaction potentials
    intPot!(poly, alpha)      # polymer self-interaction (SWEEP)
    intPot!(bulk, -1.0)       # bulk self-interaction (fixed)
    intPot!(bulk, poly, -1.0) # bulk-polymer interaction (fixed)

    # initialize lattice and species
    init!(l, [poly, bulk])

    # moving species
    movingsps = [(poly, poly.intIdx), (bulk, bulk.intIdx)]

    nmoves = [0, 0]
    niter = 0
    frame = 1

    radii = zeros(nframes)
    hams_poly = zeros(nframes)
    hams_bulk = zeros(nframes)
    accs = zeros(nframes)
    deltaH = zeros(2)
    
    # initial values
    radii[frame] = radius(poly)
    hams_poly[frame] = poly.ham
    hams_bulk[frame] = bulk.ham
    accs[frame] = 0
    
    # run simulation
    while frame <= nframes - 1
        niter += 1
        
        for (s, i) in movingsps
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH
                
                # update partner species energies
                for (p, j) in s.intSpecies
                    if j != s.intIdx
                        deltaH[j] += p.deltaHP
                    end
                end
                
                # record data at sampling intervals
                if mod(nmoves[i], step) == 0
                    frame += 1
                    
                    radii[frame] = radius(poly)
                    hams_poly[frame] = hams_poly[frame-1] + deltaH[1]
                    hams_bulk[frame] = hams_bulk[frame-1] + deltaH[2]
                    accs[frame] = nmoves[1] / niter  # track polymer acceptance
                    
                    deltaH = zeros(2)
                    
                    if frame >= nframes
                        break
                    end
                end
            end
        end
    end

    # collect results
    coords_poly = Matrix{Int32}(undef, poly.length, 4)
    coords_bulk = Matrix{Int32}(undef, bulk.length, 4)

    res["species_name"] = poly.name
    res["nframes"] = nframes
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = poly.intMat
    res["boxlength"] = poly.boxlength
    res["length"] = poly.length
    res["radius_of_gyration"] = radii
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["final_coordinates_poly"] = extractCoords!(poly, coords_poly)
    res["final_coordinates_bulk"] = extractCoords!(bulk, coords_bulk)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
