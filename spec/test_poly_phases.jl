### physical parameters

global l = Lattice(;
    boxlength = 32,
    nspecies = 2
)

# polymer
global s1 = Species(;
    length = 32,
    maxlength = 256,
    boxlength = l.boxlength,
    intPot = 1,
    chemPot = 1
)

# dilute monomers
global s2 = Species(;
    length = 256,
    maxlength = 256,
    boxlength = l.boxlength,
    intPot = 0.0,
    chemPot = 0.0
)

# initializing lattice and species
initPoly!(s1)
initDilute!(s2)
initLattice!(l, [s1, s2])

# initialize hamiltonians
initHamPoly!(s1)
hamInitial = s1.ham


### simulations parameters

niter = 1000000
step = 1000
k = 1

# alphas = ...   # parallelized parameters
# betas = ...
# nparams = length(alphas)

### functions

function target(args...)
    niter, step = args
    plotsteps = div(niter, step)
    rgs = zeros(plotsteps)

    nmoves = 0

    for i in 1:niter
        nmoves += move!(s1)
        
        if mod(i, step) == 0 && k < plotsteps
            k += 1
            rgs[k] = radius(s1)
        end
    end
end