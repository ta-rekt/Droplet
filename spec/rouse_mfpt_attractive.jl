using Droplet
import JSON
using CubicSplines: CubicSpline

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

# --- Nonlinear two-sphere potential (tabulated by the BVP solver) ----------
# Load U(d) vs surface-to-surface gap d, one curve per phi_inf value.  The
# surface-model parameters and the phi_inf sweep below are pinned to this file,
# so the linear-response and nonlinear potentials are compared at identical
# physics and each swept alpha maps one-to-one onto a stored U(d) curve.
nonlinear_path = "/home/tohme/.julia/dev/Droplet/media/figures/20260628/fisher_two_spheres_chi_2.5_h0_-0.5_R_5.0_sweep_c.json"
nl_data    = JSON.parsefile(nonlinear_path)
nl_group   = nl_data["groups"][1]
nl_d       = Vector{Float64}(nl_group["d_grid"])     # surface-to-surface gap
nl_U       = nl_group["U_full_grid"]                 # one U(d) curve per phi_inf
nl_phi_inf = Vector{Float64}(nl_group["phi_inf"])    # bulk volume fractions

# Surface adsorption model parameters (R, omega_b, chi_b pinned to the JSON).
nb              = 1.0
chi_b           = nl_data["chi_b"]
omega_b         = nl_group["omega_b"]
kappa           = 1.0
phi_inf         = 0.1
eps             = 1e-6
R               = nl_group["R"]

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.
# alphas = phi_inf, taken from the JSON so every value maps to a U(d) curve.
# beta is unused here.
alphas = nl_phi_inf
# alphas = range(0.003, stop=0.15, length=32)   # former dense linear-only sweep
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

nsteps = 10_000_000

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

    # Nonlinear two-sphere potential for this alpha: cubic-spline interpolant of
    # the tabulated U(d).  d is the surface-to-surface gap, mapped to the
    # bead-to-bead distance as d = x - R, and flat-clamped outside the data.
    nl_idx = argmin(abs.(nl_phi_inf .- alpha))
    @assert isapprox(nl_phi_inf[nl_idx], alpha; atol=1e-9) "alpha=$alpha has no matching phi_inf curve in $nonlinear_path"
    nl_U_curve = Float64.(nl_U[nl_idx])
    nl_ok      = isfinite.(nl_U_curve)        # drop the NaN at contact (d=d_grid[1])
    nl_d_ok    = nl_d[nl_ok]
    nl_spline  = CubicSpline(nl_d_ok, nl_U_curve[nl_ok])
    nl_dmin, nl_dmax = first(nl_d_ok), last(nl_d_ok)
    adsorption_potential_nonlinear(x) = nl_spline[clamp(x - p.R, nl_dmin, nl_dmax)]

    # --- Select the pair potential: keep exactly ONE definition of
    # `adsorption_potential` active (comment out the other).  Both
    # `attractive_force` and `total_potential` below read `adsorption_potential`.

    # (1) Linear-response analytic potential (surface condensation on a sphere):
    function adsorption_potential(x)
        return -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
                (p.kappa * (lam + 1/p.R) * x)
    end

    # (2) Nonlinear BVP potential (cubic spline of the JSON):
    # adsorption_potential(x) = adsorption_potential_nonlinear(x)

    # f(d) = -dU/dr, computed numerically via central differences from
    # adsorption_potential — single source of truth for the pair interaction.
    h_fd = 1e-6
    attractive_force = x ->
        -(adsorption_potential(x + h_fd) - adsorption_potential(x - h_fd)) / (2 * h_fd)

    function total_potential(x)
        return (adsorption_potential(x) 
                - 2*log(x) + 3*x^2 / (2*beadSep * kuhnLength^2))
    end


    # Sojourn (A/B dwell-time) bounds on the bead-to-bead distance, bracketing
    # the barrier ridge between the closed (contact) and open (entropic) wells.
    # A = {d ≤ a} is the closed/bound state, B = {d ≥ b} the open/unbound state;
    # the buffer (a, b) gives Schmitt-trigger hysteresis. TUNE these to the
    # barrier of `total_potential` (cf. locate_ridge in fisher_rate_rouse_mfpt.py).
    a_bound   = p.R + 0.5
    b_unbound = p.R + 2.0

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
        a               = a_bound,
        b               = b_unbound,
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
    # A/B sojourn durations (in integration steps; multiply by dt for time).
    res["sojourn_a"]        = a_bound
    res["sojourn_b"]        = b_unbound
    res["A_intervals"]      = r.A_intervals
    res["B_intervals"]      = r.B_intervals
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