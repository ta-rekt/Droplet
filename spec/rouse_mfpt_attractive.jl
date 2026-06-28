using Droplet

# Lattice species (unused -- the Rouse solver is off-lattice).
slist = []
checkpoint_path = ""

# Rouse polymer parameters
kBT = 1.0
kuhnLength = 1.0            # Kuhn length
friction = 1.0        # friction coefficient
stiffness = 2*kBT / kuhnLength^2  # bond stiffness
D = kBT / friction    # diffusion coefficient
N = 100            # chain length
firstBeadIdx = 1
beadSep = N - 1
dt = 0.001 * friction / stiffness   # 0.3 times the characteristic time

# Surface adsorption model parameters
nb              = 1.0
chi_b           = 2.5
omega_b         = -0.5
kappa           = 1.0
phi_inf         = 0.1
eps             = 1e-6
R               = 2.5

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.
# alphas: contact threshold (physical length).  beta is unused here.
alphas = range(0.003, stop=0.15, length=16)
# alphas = [0.01]
betas  = [-1.4]

nparams = length(alphas)

blurb = "Solving the Langevin equation for a Rouse polymer to
         measure the distribution of first passage times to and
         from the entropic potential well. Computing the Fisher
         information rate from these distributions.
         
         Starting now with the attractive potential from the linear
         response theory of surface condensation on a sphere.
         
         Trying to recover the Kramers exit rate from the absorption
         potential well."

nsteps = 50_000_000

function target!(res, args...)
    slist, keya, alpha, beta, blurb = args
    t0 = time()
    
    # Surface adsorption model parameters
    p = Parameters(
        nb              = nb,
        chi_b           = chi_b,
        omega_b         = omega_b,
        kappa           = kappa,
        phi_inf         = alpha,
        eps             = eps,
        R               = R,
    )
    
    # Parameters of the exponential-integral pair potential.
    d2f_phi_inf = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(d2f_phi_inf / p.kappa)
    dists = collect(range(p.R, p.L, length=500))

    # obtained from linear-response concentration profile
    function adsorption_potential(x)
        return -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
                (p.kappa * (lam + 1/p.R) * x)
    end

    # f(d) = -dU/dr, computed numerically via central differences from
    # adsorption_potential — single source of truth for the pair interaction.
    h_fd = 1e-6
    attractive_force = x ->
        -(adsorption_potential(x + h_fd) - adsorption_potential(x - h_fd)) / (2 * h_fd)

    function total_potential(x)
        return (adsorption_potential(x) 
                - 2*log(x) + 3*x^2 / (2*beadSep * kuhnLength^2))
    end


    # primary run: distance-only mode, long trajectory for contact statistics
    r = RousePolymer(;
        boltzmannConst  = 1.0, 
        temperature     = kBT,
        friction        = friction, 
        kuhnLength      = kuhnLength,
        r_min           = p.R,
        nbeads          = Int32(N), 
        nsteps          = Int32(nsteps),
        dt              = dt,
        trackedBead1    = firstBeadIdx, 
        trackedBead2    = firstBeadIdx + beadSep,
        attractiveForce = attractive_force,
        saveMode        = :distance_only, 
        seed            = keya + 42,
    )

    init!(r)
    simulate!(r)

    res["blurb"]            = blurb
    res["alpha"]            = alpha
    res["beta"]             = beta
    res["surface_params"]   = Dict("nb" => p.nb, "chi_b" => p.chi_b, 
                            "kappa" => p.kappa, "omega_b" => p.omega_b,
                            "R" => p.R, "phi_inf" => p.phi_inf)
    res["rouse_params"]     = Dict("N" => N, "kuhn_length" => kuhnLength, "kBT" => kBT,
                            "friction" => friction, "stiffness" => stiffness,
                            "D" => D, "dt" => dt, "bead_separation" => beadSep,
                            "nsteps" => nsteps)
    sample_skip = 500
    res["distance_traj"]    = r.distanceTraj[1:sample_skip:end]
    res["sample_skip"]      = sample_skip
    res["total_potential"]  = total_potential.(dists)
    res["bead_relax_time"]  = r.timeScale
    res["bead_relax_scale"] = r.posScale
    res["runtime_min"]      = (time() - t0) / 60.0

    return res
end

# deprecated: smooth crossover to flat potential near surface

# Effective pair potential of the bead-to-bead distance d:
#   U(d) = -4π R³ h0² e^{-λ(d-R)} / [κ(λ + 1/R) · d]   for d ≥ R,
# smoothly continued to a flat plateau for d < R via softplus(d-R; β)
# ≈ max(d-R, 0), which kills the 1/d and exp(-λ(d-R)) divergences at d→0.
# β sets the crossover width (~1/β); β=30 ⇒ width ≈ 0.03·b.
# β_smooth = 30.0
# function adsorption_potential(d)
#     Δ    = d - p.R
#     soft = max(Δ, 0) + log1p(exp(-β_smooth * abs(Δ))) / β_smooth   # ≈ max(Δ, 0)
#     deff = p.R + soft                                              # ≈ max(d, R)
#     return -(4π * p.R^3 * h0^2 * exp(-lam * soft)) /
#             (p.kappa * (lam + 1/p.R) * deff)
# end