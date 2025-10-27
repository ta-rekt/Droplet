### physical parameters

boxlength = 32
nspecies = 1
maxlength = 256

global l = Lattice(;
    boxlength = boxlength,
    nspecies = nspecies
)

# polymer
global s1 = Species(;
    length = boxlength,
    maxlength = maxlength,
    boxlength = boxlength,
    intPot = nothing,
    chemPot = nothing
)

# list containing all species
slist = [s1]


### simulation parameters

nframes = 2000
step = 50

alphas = [1.0]   # parallelized parameters (interaction potential)
betas = [1.0]    # chemical potential
nparams = length(alphas)


### functions

function target!(res, args...)
    # species list and sweep parameters
    slist, alpha, beta = args

    s1 = slist[1]
    s1.intPot = alpha
    s1.chemPot = beta

    # initialize lattice and species
    initPoly!(s1)
    initLattice!(l, [s1])

    # initialize hamiltonians
    initHamPoly!(s1)

    nmoves = 0
    niter = 0
    frame = 1

    radii = zeros(nframes)
    lengths = zeros(nframes)
    hams = zeros(nframes)

    hamInitial = s1.ham   # initial energy

    # run simluation
    while frame <= nframes
        
        if mod(nmoves, step) == 0

            radii[frame] = radius(s1)
            lengths[frame] = s1.length
            hams[frame] = frame == 1 ? hamInitial : hams[frame-1] + s1.deltaH

            frame += 1
        end

        niter += 1
        nmoves += move!(s1)
    end

    # collect results
    res["gyration_radius"] = radii
    res["polymer_length"] = lengths
    res["hamiltonian"] = hams

    return res
end
