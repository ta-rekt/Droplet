# test_bvp_solver.jl — solve_isolated_sphere (src/bvp_solver.jl) checked against
# the Cahn "graphical construction" for surface adsorption / prewetting, in the
# zero-curvature limit R -> infinity (here R = 1e4).
#
# THE CONSTRUCTION.  In planar geometry the isolated-surface BVP is
#
#     kappa * phi'' = dfb/dphi(phi) - mu_inf/nb,      mu_inf = nb*dfb_dphi(phi_inf)
#
# with phi -> phi_inf, phi' -> 0 in the far field.  Multiplying by phi' and
# integrating in from infinity gives the first integral
#
#     (kappa/2) * phi'(z)^2 = Delta_f(phi(z)),
#     Delta_f(phi) = fb(phi) - fb(phi_inf) - dfb_dphi(phi_inf)*(phi - phi_inf)
#
# (Delta_f is fb measured from the tangent at phi_inf, so it vanishes quadratically
# there and is >= 0 for any phi_inf outside the two-phase region).  The Robin
# adsorption BC is -kappa*phi'(0) + omega_b = 0, i.e. phi'(0) = h0/kappa, so the
# surface composition is fixed by a single scalar equation -- the horizontal line
# cutting the Delta_f curve that the notebook draws:
#
#     Delta_f(phi_s) = h0^2 / (2*kappa).                                     (*)
#
# phi_s is the SMALLEST root of (*) above phi_inf: that is the branch reached by
# integrating in from the dilute far field, and the one the solver's continuation
# tracks.  If (*) has no root at all then Delta_f < h0^2/(2 kappa) everywhere,
# phi' never vanishes, the profile runs into phi = 1 at finite z, and the surface
# is completely wet -- the solver must then either saturate or fail, never report
# a converged finite film.
#
# bc_reg BREAKS (*).  With bc_reg > 0 the surface BC carries the regularizing
# log term (see _surf_slope), phi'(0) = (h0 + bc_reg*log(phi_s/(1-phi_s)))/kappa,
# and the construction becomes
#
#     Delta_f(phi_s) = (h0 + bc_reg*log(phi_s/(1-phi_s)))^2 / (2*kappa).     (**)
#
# Below phi_s = 1/2 the log term is negative and so reinforces an attractive
# h0 < 0: a bc_reg > 0 solve sits at a systematically HIGHER phi_s than the bare
# h0^2/(2 kappa) construction predicts.  Comparing the two is a mismatch of
# boundary conditions, not a solver error.  "BcRegShiftsTheConstruction" below
# pins both forms down.
#
# FREE ENERGY.  fb here is the ASYMMETRIC (polymer-in-solvent) Flory-Huggins form
# of src/bvp_solver.jl, whose critical point is chi_c = (1 + 1/sqrt(nb))^2/2 and
# whose binodal must be found by common tangent.  It is NOT the symmetric blend
# that spec/.../param.jl's `dilute_binodal` assumes (that one is the binodal of
# an effective chi = nb*chi_b and agrees only at nb = 1), so this file computes
# its own binodal.

using Test
using Random
using Printf
using Droplet

# ---------------------------------------------------------------------------
# The graphical construction, evaluated exactly
# ---------------------------------------------------------------------------

# Bisection on a bracketed sign change; 80 halvings exhaust Float64.
function _bisect(f, a, b; n = 80)
    fa = f(a)
    for _ in 1:n
        m = 0.5 * (a + b)
        (fa * f(m) <= 0) ? (b = m) : (a = m; fa = f(m))
    end
    return 0.5 * (a + b)
end

# fb measured from its tangent at phi_inf -- the curve the construction reads off.
deltaf(phi, p::Parameters) =
    fb(phi, p) - fb(p.phi_inf, p) - dfb_dphi(p.phi_inf, p) * (phi - p.phi_inf)

chi_critical(nb) = (1 + 1 / sqrt(nb))^2 / 2

# Common-tangent binodal (phi_dilute, phi_dense) of the asymmetric fb.  dfb_dphi
# rises, falls between the two spinodal points, then rises again, so each branch
# holds exactly one root of dfb_dphi = mu; bisect on mu until the two branches
# have equal grand potential fb - mu*phi.
function binodal(chi, nb)
    p = Parameters(nb = nb, chi_b = chi, omega_b = 0.0, kappa = 1.0,
                   phi_inf = 0.1, L = 1.0, R = Inf, eps = 1e-16)
    lo_sp, hi_sp = compute_spinodal(p)
    lo_sp === nothing && return (nothing, nothing)
    dilute(mu) = _bisect(x -> dfb_dphi(x, p) - mu, 1e-300, lo_sp)
    dense(mu)  = _bisect(x -> dfb_dphi(x, p) - mu, hi_sp, 1 - 1e-16)
    gap(mu)    = (fb(dense(mu), p) - mu * dense(mu)) -
                 (fb(dilute(mu), p) - mu * dilute(mu))
    mu = _bisect(gap, dfb_dphi(hi_sp, p), dfb_dphi(lo_sp, p))
    return dilute(mu), dense(mu)
end

# Delta_f has at most two interior turning points -- the other two solutions of
# dfb_dphi(phi) = dfb_dphi(phi_inf), one in each spinodal bracket -- so it is
# monotone on three segments: rising to the (barrier) local max, falling to the
# (well) local min on the dense branch, rising again into phi = 1.
function turning_points(p::Parameters)
    lo_sp, hi_sp = compute_spinodal(p)
    slope = dfb_dphi(p.phi_inf, p)
    return (_bisect(x -> dfb_dphi(x, p) - slope, lo_sp, hi_sp),
            _bisect(x -> dfb_dphi(x, p) - slope, hi_sp, 1 - 1e-16))
end

# Smallest root of Delta_f(phi) = E(phi) above phi_inf, or NaN when there is none
# (complete wetting).  E is the level line: the constant h0^2/(2 kappa) of (*),
# or the phi_s-dependent right-hand side of (**) when bc_reg > 0.
function surface_root(p::Parameters, E)
    barrier, well = turning_points(p)
    g(x) = deltaf(x, p) - (E isa Function ? E(x) : E)
    for (a, b) in ((p.phi_inf * (1 + 1e-12), barrier), (barrier, well),
                   (well, 1 - 1e-16))
        g(a) * g(b) <= 0 && return _bisect(g, a, b)
    end
    return NaN
end

# ---------------------------------------------------------------------------
# Solver domain
# ---------------------------------------------------------------------------
# The truncation z in [0, L] must be long enough that the Dirichlet wall at z = L
# does no work: too short a box props the profile up and the solver happily
# "converges" to a finite film that the construction says cannot exist.  The
# first integral gives the exact profile width by quadrature,
# dz = dphi / sqrt(2*Delta_f(phi)/kappa) (cubic-clustered near the logarithmically
# divergent phi_inf end, cut off at phi_inf*(1+1e-3)), and the far tail decays over
# ell = sqrt(kappa/fb''(phi_inf)); L = 2*width + 20*ell covers both with margin.
# Complete-wetting cases have no phi_s, so their width runs all the way to phi = 1.
function _film_width(p::Parameters, phi_t::Real; n = 1200)
    lo = p.phi_inf * (1 + 1e-3)
    w, x_prev = 0.0, lo
    for k in 1:n
        x = lo + (phi_t - lo) * (k / n)^3
        d = max(deltaf(0.5 * (x + x_prev), p), 0.0)
        d > 0 && (w += (x - x_prev) / sqrt(2 * d / p.kappa))
        x_prev = x
    end
    return w
end

solver_domain(p::Parameters, phi_s::Real) =
    2 * _film_width(p, isnan(phi_s) ? 1 - 1e-6 : phi_s) +
    12 * sqrt(p.kappa / d2fb_dphi2(p.phi_inf, p))

# ---------------------------------------------------------------------------

@testset "bvp_solver" begin

    # -- 1. The first-integral identity on a few hand-picked points -----------
    # Weak, moderate and near-prewetting adsorption at the pipeline's own
    # chi_b = 2.5, nb = 1: phi_s must sit on Delta_f = h0^2/(2 kappa) to the
    # 1/R accuracy of the R = 1e4 sphere.
    @testset "FirstIntegral" begin
        for (h0, phi_inf) in ((-0.2, 0.05), (-0.5, 0.05), (-0.5, 0.09),
                              (-0.5, 0.12), (-0.7, 0.03))
            p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = h0, kappa = 1.0,
                           phi_inf = phi_inf, L = 1.0, R = 1e4, eps = 1e-16)
            E = h0^2 / (2 * p.kappa)
            r = surface_root(p, E)
            L = solver_domain(p, r)
            s = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                      bc_reg = 0.0)
            @test s.success
            @test isapprox(s.phi[1], r; atol = 1e-3)
            @test isapprox(deltaf(s.phi[1], p), E; atol = 1e-3)
            @test s.phi[end] ≈ phi_inf                 # Dirichlet far field
            @test all(diff(s.phi) .<= 1e-9)            # monotone decay
        end
    end

    # -- 2. bc_reg moves the construction, it does not break the solver -------
    # The regression this file was written for: spec/jobs/20260804/run1 solves
    # with bc_reg = 5e-2 and was compared against the bare h0^2/(2 kappa) line.
    # At bc_reg = 0 the solver lands on (*); at bc_reg > 0 it lands on (**), and
    # the two differ by far more than any numerical tolerance.
    @testset "BcRegShiftsTheConstruction" begin
        bc_reg = 5e-2
        for phi_inf in (0.05, 0.065, 0.0786)
            p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.5, kappa = 1.0,
                           phi_inf = phi_inf, L = 1.0, R = 1e4, eps = 1e-16)
            E = p.omega_b^2 / (2 * p.kappa)
            r = surface_root(p, E)
            L = solver_domain(p, r)

            bare = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                         bc_reg = 0.0)
            @test bare.success
            @test isapprox(bare.phi[1], r; atol = 1e-3)

            # (**): Delta_f(phi_s) = (h0 + bc_reg*log(phi_s/(1-phi_s)))^2/(2 kappa).
            reg = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 1500,
                                        bc_reg = bc_reg)
            r_reg = surface_root(p, x -> (p.omega_b +
                                          bc_reg * log(x / (1 - x)))^2 / (2 * p.kappa))
            @test reg.success
            @test isapprox(reg.phi[1], r_reg; atol = 1e-3)

            # ... and the shift is large.  These phi_s all sit below 1/2, where the
            # log term is negative and so reinforces the attractive h0: the
            # regularized surface is markedly WETTER than the bare construction.
            @test reg.phi[1] > bare.phi[1] + 0.02
        end
    end

    # -- 3. Deep robustness scan against the construction ---------------------
    # 1000 random (omega_b, chi_b, nb) x 10 far-field compositions spanning
    # 0.1-0.8 of the dilute binodal (equivalently 10 values of mu_inf, which is
    # monotone in phi_inf below the binodal).
    #
    # chi_b is drawn from [1.5, 2.5] EXCEPT that a binodal only exists above
    # chi_critical(nb); chi_critical = 2.0 at nb = 1 (and 1.05 / 0.87 at nb =
    # 5 / 10), so at nb = 1 the draw is from [2.04, 2.5] -- below that there is
    # no two-phase region and "a fraction of the binodal" is undefined.
    #
    # Every sample is checked.  Where (*) has a root the surface composition must
    # match it; where it has none the construction predicts complete wetting, and
    # the solver must not answer with a converged finite film.  Note that most of
    # the requested omega_b range is past the wetting transition: at 1000 trials
    # 4185 of the 10000 samples have a finite phi_s, 5769 are completely wet.
    #
    # ~7 min at the default 1000 trials; BVP_SCAN_TRIALS shortens it.  The seed is
    # fixed and smaller runs are a prefix of the same stream, so the tolerances
    # below (set with 2-3x headroom over the measured 1000-trial distribution)
    # stay valid at any trial count.
    @testset "ConstructionScan" begin
        ntrial = parse(Int, get(ENV, "BVP_SCAN_TRIALS", "1000"))
        nphi   = 10
        Random.seed!(20260804)

        dphi_err = Float64[]      # |phi_s(solver) - phi_s(construction)|
        energy_err = Float64[]    # |Delta_f(phi_s) - h0^2/(2 kappa)|
        n_root = n_stall = n_wet = n_wet_ok = n_false_film = 0

        for _ in 1:ntrial
            nb   = rand((1.0, 5.0, 10.0))
            cmin = max(1.5, 1.02 * chi_critical(nb))
            chi  = cmin + rand() * (2.5 - cmin)
            h0   = -2.0 * rand()
            phi_b, _ = binodal(chi, nb)

            for frac in range(0.1, 0.8; length = nphi)
                p = Parameters(nb = nb, chi_b = chi, omega_b = h0, kappa = 1.0,
                               phi_inf = frac * phi_b, L = 1.0, R = 1e4,
                               eps = 1e-16)
                E = h0^2 / (2 * p.kappa)
                r = surface_root(p, E)
                L = solver_domain(p, r)
                s = solve_isolated_sphere(with(p; L = L); L_iso = L, N = 600,
                                          maxiter = 100, bc_reg = 0.0)

                if isnan(r)                       # construction: complete wetting
                    n_wet += 1
                    s.phi[1] > 0.99 && (n_wet_ok += 1)
                    (s.success && s.phi[1] < 0.99) && (n_false_film += 1)
                elseif !s.success                 # Newton stalled; counted below
                    n_stall += 1
                else
                    n_root += 1
                    push!(dphi_err, abs(s.phi[1] - r))
                    push!(energy_err, abs(deltaf(s.phi[1], p) - E))
                end
            end
        end

        sort!(dphi_err); sort!(energy_err)
        q(v, f) = v[clamp(round(Int, f * length(v)), 1, length(v))]
        @printf("\n  ConstructionScan: %d samples | finite phi_s %d (stalled %d) | complete wetting %d\n",
                ntrial * nphi, n_root, n_stall, n_wet)
        @printf("    |dphi_s|          median %.2e  p90 %.2e  p99 %.2e  max %.2e\n",
                q(dphi_err, 0.5), q(dphi_err, 0.9), q(dphi_err, 0.99), dphi_err[end])
        @printf("    |Delta_f - h0^2/2k| median %.2e  p99 %.2e  max %.2e\n",
                q(energy_err, 0.5), q(energy_err, 0.99), energy_err[end])

        # Surface composition matches the construction on every sample (measured:
        # median 1.5e-5, p99 1.4e-3, max 3.0e-3).  The residual is the O(1/R)
        # curvature the R = 1e4 sphere still carries -- it falls off as 1/R, and at
        # R = 1e8 or in the exact planar slab the worst case drops to 3e-6.  It is
        # amplified where the h0^2/(2 kappa) line grazes the top of the Delta_f
        # barrier and the root is ill-conditioned, hence the loose ceiling against
        # a median two orders of magnitude tighter.  For scale, mistaking (**) for
        # (*) -- the bc_reg confusion of testset 2 -- costs 0.03 to 0.2.
        @test dphi_err[end] <= 1e-2
        @test q(dphi_err, 0.99) <= 3e-3
        @test q(dphi_err, 0.5) <= 1e-4

        # The same statement read off the well-conditioned axis: whatever phi_s
        # the solver returns, it sits on the Delta_f = h0^2/(2 kappa) level line.
        @test energy_err[end] <= 3e-2
        @test q(energy_err, 0.99) <= 1e-2

        # Complete wetting: the solver may saturate or stall, but it must never
        # report a converged finite adsorbed film where the construction has no
        # intersection at all.
        @test n_false_film == 0
        @test n_wet_ok >= 0.5 * n_wet

        # Newton is allowed to stall on the stiffest near-prewetting samples, but
        # must carry the overwhelming majority.
        @test n_stall <= 0.05 * (n_root + n_stall)
    end
end
