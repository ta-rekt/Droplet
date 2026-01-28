### physical parameters
# MARK: -Phys. params

l = Lattice(;
    boxlength = 42,
    nspecies = 4  # base (1), pattern1 (2), pattern2 (3), bulk (4)
)

# heterogeneous polymer with patterned regions
poly = Polymer(;
    length = 32,
    maxlength = 128,
    initShape = "line",
    nMonomerTypes = 3,  # base, pattern1, pattern2
    initPos = [5, 21, 21],
    pattern = :custom
)

# define pattern: 128 beads with two patterned regions
pattern = zeros(Int, 128)
pattern[1:28] .= 1       # base type (28 beads)
pattern[29:36] .= 2      # pattern 1 (8 beads)
pattern[37:92] .= 1      # base type (56 beads separation)
pattern[93:100] .= 3     # pattern 2 (8 beads)
pattern[101:128] .= 1    # base type (28 beads)

poly.customPattern = pattern

initMonomers!(poly)
base = poly.monomers[1]
pattern1 = poly.monomers[2]
pattern2 = poly.monomers[3]

# bulk species
bulk = FastBulkMonomer(;
    length = 4096,
    initShape = "uniform"
)

# list containing all species
slist = [base, pattern1, pattern2, bulk]


### simulation parameters
# MARK: -Sim. params

t_max = 20    # runtime in minutes (initial test)
step = 10000
t_burn_in = 0.8  # 12 seconds

vals = range(-2.0, stop=0.0, step=0.25)   # parallelized parameters (inter-pattern interaction)
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
    t_current = 0

    # species list and sweep parameters
    slist, alpha, beta = args

    base, pattern1, pattern2, bulk = slist
    
    # set interaction potentials (4x4 interaction matrix)
    # base (type 1)
    intPot!(base, 0.5)               # base self-interaction
    intPot!(base, pattern1, 0.5)     # base-pattern1
    intPot!(base, pattern2, 0.5)     # base-pattern2
    intPot!(base, bulk, 0.0)        # base-bulk (weak)
    
    # pattern1 (type 2)
    intPot!(pattern1, 0.5)           # pattern1 self-interaction
    intPot!(pattern1, pattern2, alpha)  # pattern1-pattern2 (SWEEP)
    intPot!(pattern1, bulk, 0.0)    # pattern1-bulk (strong)
    
    # pattern2 (type 3)
    intPot!(pattern2, 0.5)           # pattern2 self-interaction
    intPot!(pattern2, bulk, 0.0)    # pattern2-bulk (strong)
    
    # bulk (type 4)
    intPot!(bulk, -1.5)              # bulk self-interaction

    # initialize lattice and species
    init!(l, [base, pattern1, pattern2, bulk])

    # moving species
    movingsps = [(base, base.intIdx), (pattern1, pattern1.intIdx), 
                 (pattern2, pattern2.intIdx), (bulk, bulk.intIdx)]

    nmoves = 0
    niter = 0
    frame = 1

    radii = []
    lengths = []
    pattern_dists = []
    contact_events = []
    first_contact_frame = -1
    hams_poly = []
    hams_bulk = []
    accs = []
    # deltaH needs to be sized for the max intIdx, not l.nspecies
    max_intIdx = maximum([base.intIdx, pattern1.intIdx, pattern2.intIdx, bulk.intIdx])
    deltaH = zeros(max_intIdx)
    
    # helper function to calculate pattern center positions
    function get_pattern_centers()
        # pattern1: beads 29-36 (bead positions in linked list)
        # pattern2: beads 93-100
        
        if poly.length < 100
            # patterns not yet grown
            return nothing, nothing, Inf
        end
        
        # calculate center of mass for each pattern
        center1 = [0.0, 0.0, 0.0]
        center2 = [0.0, 0.0, 0.0]
        
        # iterate through polymer linked list to find pattern beads
        key = poly.head
        bead_num = 1
        
        for i in 1:poly.length
            pos = poly.pos[key][1:3]
            
            # check if this bead is in pattern1 (beads 29-36)
            if bead_num >= 29 && bead_num <= 36
                center1 .+= pos
            end
            
            # check if this bead is in pattern2 (beads 93-100)
            if bead_num >= 93 && bead_num <= 100
                center2 .+= pos
            end
            
            key = poly.pos[key][5]  # next bead
            bead_num += 1
        end
        
        center1 ./= 8
        center2 ./= 8
        
        # euclidean distance
        dist = sqrt(sum((center1 .- center2).^2))
        
        return center1, center2, dist
    end
    
    # helper function to check if patterns are in contact
    function patterns_in_contact()
        if poly.length < 100
            return false
        end
        
        # collect positions of pattern1 and pattern2 beads
        pattern1_positions = []
        pattern2_positions = []
        
        key = poly.head
        bead_num = 1
        
        for i in 1:poly.length
            pos = poly.pos[key][1:3]
            
            if bead_num >= 29 && bead_num <= 36
                push!(pattern1_positions, copy(pos))
            end
            
            if bead_num >= 93 && bead_num <= 100
                push!(pattern2_positions, copy(pos))
            end
            
            key = poly.pos[key][5]
            bead_num += 1
        end
        
        # check all pairs of beads between pattern1 and pattern2
        for pos1 in pattern1_positions
            for pos2 in pattern2_positions
                # check if adjacent (within 1 lattice site in each dimension, periodic BC)
                dx = abs(pos1[1] - pos2[1])
                dy = abs(pos1[2] - pos2[2])
                dz = abs(pos1[3] - pos2[3])
                
                # account for periodic boundaries
                dx = min(dx, l.boxlength - dx)
                dy = min(dy, l.boxlength - dy)
                dz = min(dz, l.boxlength - dz)
                
                if dx <= 1 && dy <= 1 && dz <= 1 && (dx + dy + dz) > 0
                    return true
                end
            end
        end
        
        return false
    end
    
    # initial values
    push!(radii, radius(poly))
    push!(lengths, poly.length)
    _, _, initial_dist = get_pattern_centers()
    push!(pattern_dists, initial_dist)
    push!(contact_events, patterns_in_contact() ? 1 : 0)
    push!(hams_poly, poly.ham)  # polymer hamiltonian is stored in poly, not individual monomers
    push!(hams_bulk, bulk.ham)
    push!(accs, 0)
    
    # run simulation
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        niter += 1

        for (s, i) in movingsps
            if move!(s)
                nmoves += 1

                # record data
                if mod(nmoves, step) == 0
                    deltaH[i] += getDeltaHam(s)

                    # update deltaH for the partners
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += getDeltaHP(p)
                        end
                    end

                    ham!(s)

                    # compute polymer hamiltonian (stored in poly, not individual monomers)
                    if s.intIdx in [base.intIdx, pattern1.intIdx, pattern2.intIdx]
                        push!(hams_poly, poly.ham)

                    elseif s.intIdx == bulk.intIdx
                        push!(hams_bulk, getHam(s))
                    end

                    push!(radii, radius(poly))
                    push!(lengths, poly.length)
                    
                    # calculate pattern distance
                    _, _, dist = get_pattern_centers()
                    push!(pattern_dists, dist)
                    
                    # check contact
                    in_contact = patterns_in_contact()
                    push!(contact_events, in_contact ? 1 : 0)
                    
                    # record first contact frame
                    if in_contact && first_contact_frame == -1
                        first_contact_frame = frame
                    end

                    # acceptance rate
                    push!(accs, nmoves / (niter * length(movingsps)))

                    frame += 1
                end
            end
        end
    end

    # save results
    res["species_name"] = "poly"
    res["nframes"] = length(radii)
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = [
        [0.5, 0.5, 0.5, -0.1],      # base
        [0.5, 0.5, alpha, -1.0],    # pattern1
        [0.5, alpha, 0.5, -1.0],    # pattern2
        [-0.1, -1.0, -1.0, -1.5]    # bulk
    ]
    res["boxlength"] = l.boxlength
    res["length_poly"] = lengths
    res["length_bulk"] = bulk.length
    res["radius_of_gyration"] = radii
    res["pattern_distances"] = pattern_dists
    res["contact_events"] = contact_events
    res["first_contact_frame"] = first_contact_frame
    res["hamiltonian_poly"] = hams_poly
    res["hamiltonian_bulk"] = hams_bulk
    res["acceptance_rate"] = accs
    
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
