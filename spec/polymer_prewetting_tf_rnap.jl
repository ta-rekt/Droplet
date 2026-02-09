### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 2
)

poly = Polymer(;
    length = 10,
    maxlength = 100,
    initShape = "line",
    nMonomerTypes = 1,  # base, base
    initPos = [14, 21, 21],
    pattern = :uniform
)

initMonomers!(poly, defaultChemPot = 1.0)

# bulk species
bulk = FastBulkMonomer(;
    length = 1,
    initShape = "uniform"
)

# list containing all species
slist = [poly, bulk]


### simulation parameters
# MARK: -Sim. params

blurb = """
Prewetting phase diagram with a crumpled polymer as the nucleating site.
"""

t_max = 60    # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

vals = range(0.02, stop=0.03, length=16)  # parallelized parameters (concentration)
alphas = zeros(16)
nrepeats = 1

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

alphas = [0.024]    # concentration
betas = [0.66]    # temperature
nparams = length(alphas)


### functions
# MARK: -Functions

function target!(res, args...)
    t_start = time_ns()  # start timer (monotonic)
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta, blurb = args

    poly, bulk = slist

    base = poly.monomers[1]
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base, 0.5)              # base self-interaction
    intPot!(base, bulk, 0.0)        # base-bulk (weak)
    
    # bulk (type 4)
    intPot!(bulk, -1.0 / beta)              # bulk self-interaction

    bulk.length = Int32(round(l.boxlength^3 * alpha))

    updateMonomers!(poly, base.intIdx)

    # initialize lattice and species
    init!(l, [poly, bulk])

    # Initialize deltaHP for both polymers (needed for multi-polymer interactions)
    poly.deltaHP = 0.0

    # moving species
    movingsps = [(poly, bulk.intIdx + 1), 
                 (bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1
    minsize = 20

    radii = []
    lengths = []
    hams_poly = []
    hams_bulk = []
    max_droplet_size = []
    accs = []
    contact_probs = []
    contact_events = []
    
    # initial values
    push!(radii, radius(poly))
    push!(lengths, poly.length)
    push!(hams_poly, poly.ham)
    push!(hams_bulk, bulk.ham)
    push!(accs, 0)
    push!(contact_probs, 0.0)
    push!(contact_events, 0)

    growPolymer!(poly)

    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0
                    ham!(s)

                    ccs = cc(bulk)
                    droplets = filter(x -> length(x) >= minsize, ccs)
                    droplet_sizes = map(x -> length(x), droplets)

                    # compute polymer hamiltonian (stored in poly, not monomers)
                    if i == bulk.intIdx + 1
                        push!(hams_poly, poly.ham)

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    # compute contact probability between base and bulk
                    pattern_length = 10
                    contacts = 0
                    
                    # walk backwards from tail to count contacts
                    current_key = poly.head
                    count = 0
                    while count < pattern_length
                        # check if this bead is base type
                        if poly.pos[current_key][6] == base.intIdx
                            # count nearest neighbors with bulk
                            contacts += nnp(base, bulk, poly.pos[current_key][1:3]) / 4
                            # move to previous bead
                            current_key = poly.pos[current_key][5]
                            count += 1
                        else
                            break  # reached end of pattern region
                        end
                    end
                    
                    push!(contact_events, contacts > 0 ? 1 : 0)
                    push!(contact_probs, contacts / pattern_length)
                    push!(max_droplet_size, maximum(droplet_sizes, init=1))
                    push!(radii, radius(poly))
                    push!(lengths, poly.length)

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

    # save results
    res["blurb"] = blurb
    res["species_name"] = "poly"
    res["nframes"] = length(radii)
    res["step"] = step
    res["temperature"] = beta
    res["concentration"] = alpha
    res["interaction_matrix"] = bulk.intMat
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_bulk"] = bulk.length
    res["radius_of_gyration"] = radii
    res["max_droplet_size"] = max_droplet_size
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    res["contact_probability"] = contact_probs
    
    # extract polymer coordinates into a matrix
    poly_coords = zeros(Int32, poly.length, 4)
    extractCoords!(poly, poly_coords)
    res["final_coordinates_poly"] = poly_coords
    
    # extract bulk coordinates into a matrix
    bulk_coords = zeros(Int32, bulk.length, 4)
    extractCoords!(bulk, bulk_coords)
    res["final_coordinates_bulk"] = bulk_coords
    
    res["runtime_min"] = t_current
    
    return res
end
