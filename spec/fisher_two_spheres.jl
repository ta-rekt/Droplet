# Fisher information rate from the full 2D two-sphere BVP, with BOTH sides of
# the capillary-bridging bifurcation.
#
# Parallelised axis: alphas = phi_inf grid.  For each phi_inf the target solves
# the axisymmetric two-sphere BVP at a grid of separations d on TWO branches:
#   * unbridged (dilute gap)  -- marched inward from large d;
#   * bridged   (dense gap)   -- marched outward from contact, seeded by a dense
#                                blob (capillary condensation).
# The equilibrium pair potential is U_full(d) = min(F_branch)(d) - F(d_ref), and
# the surface volume fraction Phi_surf(d) is read off the equilibrium branch.
# These are integrated against the polymer P_eq(d) to obtain phi_probe; the
# Fisher rate is assembled across phi_inf in the dashboard.
#
# Bridging only exists when the bulk is near saturation, so phi_inf is swept just
# BELOW the dilute binodal (set chi_b = 2.5; binodal ~ 0.145).  At the spec's old
# chi_b = 4 / phi_inf = 0.01 the bulk is far undersaturated and no bridged branch
# exists -- see scouting notes.

# Driver: sim/sweep.jl.  Pipe builder: runJL.

using Droplet

# Lattice species (unused -- the BVP doesn't simulate particles).
slist = []
checkpoint_path = ""

# P_strong3() = Parameters(chi_b = 4.0, omega_b = -1.35, kappa = 1.0,
#                         phi_inf = 0.01, R = 10.0, eps = 1e-6)    

# (name = "strong3", p = P_strong3(), bc_reg = 0.05, d = 10.0,
#     Nrho = 467, dz_target = 0.03, rho_buf = 4.0, z_buf = 3.0,
#     N1d = 1400, L_iso = 12.0, warm = true,  tol_prof = 0.006, tol_iso = 0.01),

# Fixed physical and numerical configuration.
const R_FIXED       = 5.0
const NB            = 1.0
const CHI_B         = 2.5            # > 2 so bulk coexistence (bridging) exists
const KAPPA         = 1.0
const L_PARAM       = 20.0
const EPS_PARAM     = 1e-6

const D_MIN_FACTOR  = 0.0            # include contact (d = 0) exactly
const D_MAX_FACTOR  = 0.5
const N_D           = 30
const D_REF_FACTOR  = 8.0

const NRHO          = 200
const DZ_TARGET     = 0.05
const RHO_MAX       = 15.0
const Z_BUF         = 3.0
const F_TOL         = 1e-6
const MAXITER       = 200
const BRIDGE_PHI    = 0.5            # gap-centre phi above this => "bridged"
    
# Polymer (matches fisher_rate_from_profile)
const S_POLYMER     = 150.0
const A_BOND        = 1.0

# Dilute binodal phi_b(chi): small root of the symmetric-FH coexistence
# condition log(phi/(1-phi)) + nb*chi*(1-2phi) = 0.  Bridging needs phi_inf just
# below this; the dense coexistence value is 1 - phi_b (used to seed the bridge).
function _dilute_binodal(chi::Real; nb::Real = 1.0)
    g(phi) = log(phi / (1.0 - phi)) + nb * chi * (1.0 - 2.0 * phi)
    lo, hi = 1e-8, 0.5 - 1e-9
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if g(mid) < 0.0; lo = mid; else; hi = mid; end
    end
    return 0.5 * (lo + hi)
end
const BINODAL = _dilute_binodal(CHI_B; nb = NB)

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.  phi_inf is swept
# as a fraction of the binodal (near-saturation), where bridging is active.
betas = BINODAL .* collect(range(0.05, 0.90; length = 32))   # phi_inf
# betas = [BINODAL * 0.5]
alphas  = [-0.4]                                                  # omega_b

nparams = length(alphas)

blurb = "Fisher information rate from the full 2D two-sphere BVP, capturing " *
        "both sides of the capillary-bridging bifurcation.  Sweeps phi_inf " *
        "(alphas, as a fraction of the dilute binodal) at fixed (omega_b, R, " *
        "chi_b=2.5); for each phi_inf solves the unbridged (inward) and " *
        "bridged (outward, dense-seeded) branches on a d-grid, takes the " *
        "equilibrium U_full(d)=min(F)-F_ref and Phi_surf(d) on that branch, " *
        "and integrates the polymer P_eq(d) to get phi_probe."

# Fixed, anchored grid.  (1) build_geometry spans z over length 4R + d + 2*z_buf;
# (2) pin dz EXACTLY to DZ_TARGET and place z_buf on the grid so the symmetry
# plane (sphere 1 at z=0) is a grid line; snapping d to a multiple of dz puts
# BOTH sphere surfaces on grid lines for every d, so the film free energy
# discretises identically and cancels in U = F(d) - F_ref.
const Z_BUF_GRID = ceil(Int, (R_FIXED + Z_BUF) / DZ_TARGET) * DZ_TARGET - R_FIXED
_snap(d::Real)   = round(d / DZ_TARGET) * DZ_TARGET
_nz_for(d::Real) = Int(round((4.0 * R_FIXED + _snap(d) + 2.0 * Z_BUF_GRID) / DZ_TARGET)) + 1

# Warm-start guess for `geom` from a previous solution (phi_prev on geom_prev):
# sample it at geom's unknown cells, folding sphere-2's half onto the fixed
# sphere-1 half (sphere 1 stays at z=0 for all d; sphere 2 is its mirror image).
function _warm_x0(geom, phi_prev, geom_prev)
    z_mid = 0.5 * (geom.z_c1 + geom.z_c2)
    rho_q = Float64[]
    z_q   = Float64[]
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        geom.unk_mask[i, j] || continue
        push!(rho_q, geom.rho_grid[i])
        zz = geom.z_grid[j]
        push!(z_q, zz > z_mid ? 2.0 * z_mid - zz : zz)
    end
    return clamp.(interp_phi(phi_prev, geom_prev, rho_q, z_q),
                  geom.p.eps, 1.0 - geom.p.eps)
end

# Bridged seed: dense phase within `reach` of either sphere surface (fills the
# neck), bulk phi_inf elsewhere.  Bootstraps the capillary-bridge branch.
function _dense_seed(geom, phi_dense; reach = 2.0)
    R = geom.p.R
    x = Float64[]
    @inbounds for j in 1:geom.Nz, i in 1:geom.Nrho
        geom.unk_mask[i, j] || continue
        rho = geom.rho_grid[i]; z = geom.z_grid[j]
        s = min(sqrt(rho^2 + (z - geom.z_c1)^2),
                sqrt(rho^2 + (z - geom.z_c2)^2)) - R
        push!(x, s < reach ? phi_dense : geom.p.phi_inf)
    end
    return clamp.(x, geom.p.eps, 1.0 - geom.p.eps)
end

# Gap-centre phi (on axis, at the midplane) -- the bridged/unbridged label.
_gap_center(geom, phi) =
    phi[1, argmin(abs.(geom.z_grid .- 0.5 * (geom.z_c1 + geom.z_c2)))]

function target!(res, args...)
    slist, keya, alpha, beta, blurb = args

    t0 = time()

    omega_b = alpha
    phi_inf = beta
    R       = R_FIXED

    p = Parameters(nb              = NB,
                   chi_b           = CHI_B,
                   kappa           = KAPPA,
                   L               = L_PARAM,
                   eps             = EPS_PARAM,
                   conserve_volume = false,
                   parallel_plates = false,
                   omega_b         = omega_b,
                   R               = R,
                   phi_inf         = phi_inf)

    d_min  = D_MIN_FACTOR * R
    d_max  = D_MAX_FACTOR * R
    d_grid = _snap.(collect(range(d_min, d_max; length = N_D)))
    d_ref  = _snap(D_REF_FACTOR * R)
    phi_dense = 1.0 - BINODAL

    geom_at(dd) = build_geometry(p, dd; rho_max = RHO_MAX, z_buf = Z_BUF_GRID,
                                 Nrho = NRHO, Nz = _nz_for(dd))
    _solve(dd; x0 = nothing, omega_cont = true) =
        solve_two_spheres(p, dd; rho_max = RHO_MAX, z_buf = Z_BUF_GRID,
                          Nrho = NRHO, Nz = _nz_for(dd),
                          f_tol = F_TOL, maxiter = MAXITER,
                          x0_init = x0, omega_continuation = omega_cont,
                          verbose = false, bc_reg = 5e-2)

    # Reference: spheres far apart on the dilute branch.
    res_ref = _solve(d_ref)
    F_ref = res_ref[:success] ? excess_free_energy_two_spheres(res_ref) : NaN

    res_unbr = Vector{Any}(nothing, N_D)   # dilute branch result per d
    res_brid = Vector{Any}(nothing, N_D)   # bridged branch result per d
    F_unbr   = fill(NaN, N_D)
    F_brid   = fill(NaN, N_D)

    # --- Unbridged (dilute) branch: inward march from large d ---
    prev_phi  = res_ref[:success] ? res_ref[:phi]  : nothing
    prev_geom = res_ref[:success] ? res_ref[:geom] : nothing
    for k in sortperm(d_grid; rev = true)
        dd = d_grid[k]
        x0 = prev_phi === nothing ? nothing :
             _warm_x0(geom_at(dd), prev_phi, prev_geom)
        rk = _solve(dd; x0 = x0)
        (rk[:success] && _gap_center(rk[:geom], rk[:phi]) < BRIDGE_PHI) || continue
        res_unbr[k] = rk
        F_unbr[k]   = excess_free_energy_two_spheres(rk)
        prev_phi    = rk[:phi]
        prev_geom   = rk[:geom]
    end

    # --- Bridged (dense gap) branch: outward march from contact ---
    prev_phi = nothing; prev_geom = nothing
    for k in sortperm(d_grid)
        dd = d_grid[k]
        geom_k = geom_at(dd)
        guesses = prev_phi === nothing ?
            (_dense_seed(geom_k, phi_dense),) :
            (_warm_x0(geom_k, prev_phi, prev_geom), _dense_seed(geom_k, phi_dense))
        for x0 in guesses
            rk = _solve(dd; x0 = x0, omega_cont = false)
            if rk[:success] && _gap_center(rk[:geom], rk[:phi]) > BRIDGE_PHI
                res_brid[k] = rk
                F_brid[k]   = excess_free_energy_two_spheres(rk)
                prev_phi    = rk[:phi]
                prev_geom   = rk[:geom]
                break
            end
        end
    end

    # --- Backfill: inward march over any d the outward march left unsolved.
    # Cold dense seeds are fragile at isolated d; warm-starting from the
    # nearest converged bridged solution at larger d usually converges. ---
    prev_phi = nothing; prev_geom = nothing
    for k in sortperm(d_grid; rev = true)
        if res_brid[k] !== nothing
            prev_phi  = res_brid[k][:phi]
            prev_geom = res_brid[k][:geom]
            continue
        end
        prev_phi === nothing && continue
        dd = d_grid[k]
        rk = _solve(dd; x0 = _warm_x0(geom_at(dd), prev_phi, prev_geom),
                    omega_cont = false)
        if rk[:success] && _gap_center(rk[:geom], rk[:phi]) > BRIDGE_PHI
            res_brid[k] = rk
            F_brid[k]   = excess_free_energy_two_spheres(rk)
            prev_phi    = rk[:phi]
            prev_geom   = rk[:geom]
        end
    end

    # --- Equilibrium = lower-free-energy branch at each d ---
    U_unbr   = fill(NaN, N_D)
    U_brid   = fill(NaN, N_D)
    U_full   = fill(NaN, N_D)      # equilibrium pair potential
    F_arr    = fill(NaN, N_D)
    Phi_surf = fill(NaN, N_D)
    branch   = zeros(Int, N_D)     # 0 none, 1 unbridged, 2 bridged
    south    = [Float64(π)]

    # Per-d full 2D phi(rho,z) of the equilibrium branch.  Stored as
    # Vector{Vector{Float64}}, outer index rho, inner z, so on the Python side
    # `np.array(profile_phi[k])` is shape (Nrho, Nz).
    profile_d   = Float64[]
    profile_rho = Vector{Vector{Float64}}()
    profile_z   = Vector{Vector{Float64}}()
    profile_phi = Vector{Vector{Vector{Float64}}}()

    for k in 1:N_D
        fu = F_unbr[k]; fb = F_brid[k]
        if isfinite(F_ref)
            isnan(fu) || (U_unbr[k] = fu - F_ref)
            isnan(fb) || (U_brid[k] = fb - F_ref)
        end
        chosen = nothing; F_eq = NaN
        if !isnan(fu) && (isnan(fb) || fu <= fb)
            chosen = res_unbr[k]; F_eq = fu; branch[k] = 1
        elseif !isnan(fb)
            chosen = res_brid[k]; F_eq = fb; branch[k] = 2
        end
        chosen === nothing && continue
        F_arr[k]    = F_eq
        U_full[k]   = isfinite(F_ref) ? F_eq - F_ref : NaN
        Phi_surf[k] = phi_on_sphere(chosen[:phi], chosen[:geom],
                                    chosen[:geom].z_c1, south)[1]
        push!(profile_d,   d_grid[k])
        push!(profile_rho, copy(chosen[:geom].rho_grid))
        push!(profile_z,   copy(chosen[:geom].z_grid))
        push!(profile_phi, [chosen[:phi][i, :] for i in 1:size(chosen[:phi], 1)])
    end

    ok = isfinite.(U_full) .& isfinite.(Phi_surf)
    n_ok = count(ok)
    phi_probe = NaN
    success   = false
    if n_ok >= 2
        d_ok   = d_grid[ok]
        U_ok   = U_full[ok]
        Phi_ok = Phi_surf[ok]
        U_ent  = 3.0 .* (d_ok .+ R) .^ 2 ./ (2.0 * S_POLYMER * A_BOND^2)
        P_eq_un = (d_ok .+ R) .^ 2 .* exp.(-U_ent .- U_ok)
        denom   = sum(P_eq_un)
        if denom > 0.0 && isfinite(denom)
            P_eq    = P_eq_un ./ denom
            phi_probe = sum(P_eq .* Phi_ok)
            success   = isfinite(F_ref) && all(ok)
        end
    end

    res["blurb"]         = blurb
    res["phi_inf"]       = phi_inf
    res["omega_b"]       = omega_b
    res["R"]             = R
    res["chi_b"]         = CHI_B
    res["binodal"]       = BINODAL
    res["s"]             = S_POLYMER
    res["a"]             = A_BOND
    res["d_grid"]        = d_grid
    res["d_ref"]         = d_ref
    res["F_ref"]         = F_ref
    res["F"]             = F_arr
    res["U_full"]        = U_full          # equilibrium (min-F) branch
    res["U_unbridged"]   = U_unbr
    res["U_bridged"]     = U_brid
    res["branch"]        = branch          # 1 unbridged, 2 bridged, 0 none
    res["Phi_surf"]      = Phi_surf
    res["phi_probe"]     = phi_probe
    res["success"]       = success
    res["n_d_failed"]    = N_D - n_ok
    res["runtime_min"]   = (time() - t0) / 60.0
    res["profile_d"]     = profile_d
    res["profile_rho"]   = profile_rho
    res["profile_z"]     = profile_z
    res["profile_phi"]   = profile_phi
    return res
end
