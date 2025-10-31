### physical parameters

boxlength = 32
nspecies = 1
maxlength = 512

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

nframes = 500000
step = 10000

alphas = range(-1, stop=1, step=0.25)   # parallelized parameters (interaction potential)
betas = [1.0]    # chemical potential
nparams = length(alphas)


### functions

function target!(out, args...)
    t_start = time_ns()  # start timer (monotonic)

    # species list and sweep parameters
    slist, ka, alpha, kb, beta = args
    s1 = slist[1]
    s1.intPot = alpha
    s1.chemPot = beta

    # initialize lattice and species
    initPoly!(s1)
    initLattice!(l, [s1])

    # initialize hamiltonians
    initHamPoly!(s1)
    hamInitial = s1.ham   # initial energy

    nmoves = 0
    niter = 0
    frame = 1
    ham = 0

    # ensure outer-key exists and is a dict
    if !haskey(out, ka) || !(out[ka] isa AbstractDict)
        out[ka] = Dict{Any,Any}()
    end

    # ensure inner-key exists and is a dict
    if !haskey(out[ka], kb) || !(out[ka][kb] isa AbstractDict)
        out[ka][kb] = Dict{Any,Any}()
    end
    
    # collect output metadata
    coords_frame = Matrix{Int32}(undef, s1.maxlength, 3)
    out[ka][kb]["species_name"] = s1.name
    out[ka][kb]["nframes"] = nframes
    out[ka][kb]["step"] = step
    out[ka][kb]["alpha"] = alpha
    out[ka][kb]["beta"] = beta
    out[ka][kb]["maxlength"] = s1.maxlength
    out[ka][kb]["boxlength"] = s1.boxlength
    
    # initialize arrays for simulation data
    out[ka][kb]["gyration_radius"] = zeros(nframes)
    out[ka][kb]["polymer_length"] = zeros(nframes)
    out[ka][kb]["hamiltonian"] = zeros(nframes)
    out[ka][kb]["acceptance_rate"] = zeros(nframes)
    out[ka][kb]["coordinates"] = zeros(nframes, s1.maxlength, 3)
    out[ka][kb]["timestamp"] = zeros(nframes)

    # run simluation
    while frame <= nframes

        if mod(nmoves, step) == 0

            # update Hamiltonian
            ham = frame == 1 ? hamInitial : ham + s1.deltaH

            out[ka][kb]["gyration_radius"][frame] = radius(s1)
            out[ka][kb]["polymer_length"][frame] = s1.length
            out[ka][kb]["hamiltonian"][frame] = ham
            out[ka][kb]["acceptance_rate"][frame] = niter == 0 ? 0 : nmoves / niter
            out[ka][kb]["coordinates"][frame] = extractCoords!(s1, coords_frame)
            out[ka][kb]["timestamp"][frame] = (time_ns() - t_start) / 1e9

            frame += 1

            # write data
            open(file, "w") do f
                JSON.print(f, out, 4)
            end
        end
        
        niter += 1
        nmoves += move!(s1)
    end

    return out
end
