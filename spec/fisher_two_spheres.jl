# Fisher information rate from the full 2D two-sphere BVP.
#
# Parallelised axis: alphas = phi_inf grid.  For each phi_inf the target
# solves the axisymmetric two-sphere BVP at a grid of separations d, records
# the pair potential U_full(d) = F(d) - F(d_ref) and the surface volume
# fraction Phi_surf(d) at sphere 1's south pole (outward-facing, on the
# symmetry axis, toward the z end-cap), and integrates them against the
# polymer P_eq(d) to obtain phi_probe.  The Fisher rate is assembled across
# phi_inf in the dashboard (it needs the full array).
#
# Driver: sim/sweep.jl.  Pipe builder: runJL.
#
# Julia port of spec/fisher_two_spheres.py.

using Droplet

# Lattice species (unused -- the BVP doesn't simulate particles).
slist = []
checkpoint_path = ""

# Fixed physical and numerical configuration (matches the Python spec).
const R_FIXED       = 5.0
const NB            = 1.0
const CHI_B         = 4.0
const KAPPA         = 1.0
const L_PARAM       = 20.0
const EPS_PARAM     = 1e-6

const D_MIN_FACTOR  = 0.05
const D_MAX_FACTOR  = 3.0
const N_D           = 20
const D_REF_FACTOR  = 8.0

const NRHO          = 100
const DZ_TARGET     = 0.05
const RHO_MAX       = 10.0
const Z_BUF         = 4.0
const F_TOL         = 1e-6
const MAXITER       = 200

# Polymer (matches fisher_rate_from_profile)
const S_POLYMER     = 150.0
const A_BOND        = 1.0

# Sweep axes.  alphas is sharded across CPUs by sim/sweep.jl.
alphas = collect(range(0.001, 0.01; length = 20))   # phi_inf
betas  = [-1.5]                                     # omega_b

nparams = length(alphas)

blurb = "Fisher information rate from the full 2D two-sphere BVP.  " *
        "Sweeps phi_inf (alphas) at fixed (omega_b, R); for each phi_inf " *
        "solves the 2D axisymmetric BVP on a d-grid, records U_full(d) " *
        "and Phi_surf(d) at sphere 1's outward-facing pole, and " *
        "integrates the polymer P_eq(d) to get phi_probe.  Fisher rate " *
        "is computed downstream by the dashboard."

_nz_for(d::Real) = Int(floor((2.0 * R_FIXED + d + 2.0 * Z_BUF) / DZ_TARGET)) + 1

function target!(res, args...)
    slist, keya, alpha, beta, blurb = args

    t0 = time()

    omega_b = beta
    phi_inf = alpha
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
    d_grid = collect(range(d_min, d_max; length = N_D))
    d_ref  = D_REF_FACTOR * R

    function _solve(dd)
        return solve_two_spheres(p, dd;
                                 rho_max = RHO_MAX, z_buf = Z_BUF,
                                 Nrho    = NRHO,   Nz    = _nz_for(dd),
                                 f_tol   = F_TOL,  maxiter = MAXITER,
                                 verbose = false)
    end

    res_ref = _solve(d_ref)
    F_ref = res_ref[:success] ?
        excess_free_energy_two_spheres(res_ref) : NaN

    U_full   = fill(NaN, N_D)
    Phi_surf = fill(NaN, N_D)
    F_arr    = fill(NaN, N_D)
    south    = [Float64(π)]

    for k in 1:N_D
        dd = d_grid[k]
        rk = _solve(dd)
        rk[:success] || continue
        Fk = excess_free_energy_two_spheres(rk)
        F_arr[k]    = Fk
        U_full[k]   = isfinite(F_ref) ? Fk - F_ref : NaN
        Phi_surf[k] = phi_on_sphere(rk[:phi], rk[:geom], rk[:geom].z_c1, south)[1]
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
            success   = all(ok) && isfinite(F_ref)
        end
    end

    res["blurb"]       = blurb
    res["phi_inf"]     = phi_inf
    res["omega_b"]     = omega_b
    res["R"]           = R
    res["chi_b"]       = CHI_B
    res["s"]           = S_POLYMER
    res["a"]           = A_BOND
    res["d_grid"]      = d_grid
    res["d_ref"]       = d_ref
    res["F_ref"]       = F_ref
    res["F"]           = F_arr
    res["U_full"]      = U_full
    res["Phi_surf"]    = Phi_surf
    res["phi_probe"]   = phi_probe
    res["success"]     = success
    res["n_d_failed"]  = N_D - n_ok
    res["runtime_min"] = (time() - t0) / 60.0
    return res
end
