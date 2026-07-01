# Flory-Huggins bulk thermodynamics + 1D isolated-sphere BVP solver.
#
# Mirrors viz_util/bvp_solver.py: provides the Parameters dataclass and the
# scalar thermodynamic functions (fb, dfb_dphi, d2fb_dphi2, psi_linear,
# compute_mu_inf, compute_spinodal) that the 2D two-sphere code calls into.
#
# Only the minimum 1D solver needed by the two-sphere test suite is ported
# (the test's "single-sphere limit" comparison).  The full robust/shooting/
# continuation pipeline in the Python file is not reproduced.

is_spherical(p::Parameters) = isfinite(p.R)

function with(p::Parameters;
              nb=p.nb, chi_b=p.chi_b, omega_b=p.omega_b, kappa=p.kappa,
              phi_inf=p.phi_inf, L=p.L, R=p.R, eps=p.eps,
              conserve_volume=p.conserve_volume,
              parallel_plates=p.parallel_plates)
    Parameters(; nb=nb, chi_b=chi_b, omega_b=omega_b, kappa=kappa,
                 phi_inf=phi_inf, L=L, R=R, eps=eps,
                 conserve_volume=conserve_volume,
                 parallel_plates=parallel_plates)
end

@inline clip01(x::Real, eps_::Real=1e-10) = clamp(x, eps_, 1.0 - eps_)

@inline function regularized_log(x::Real, eps_::Real=1e-10)
    x > eps_ ? log(x) : log(eps_) + (x - eps_) / eps_
end

@inline function fb(phi::Real, p::Parameters)
    phi_c = clip01(phi, p.eps)
    (1.0 / p.nb) * (phi_c * regularized_log(phi_c, p.eps) +
                    (1 - phi_c) * regularized_log(1 - phi_c, p.eps)) +
        p.chi_b * phi_c * (1 - phi_c)
end

@inline function dfb_dphi(phi::Real, p::Parameters)
    phi_c = clip01(phi, p.eps)
    (1.0 / p.nb) * (regularized_log(phi_c, p.eps) - regularized_log(1 - phi_c, p.eps)) +
        p.chi_b * (1.0 - 2.0 * phi_c)
end

@inline function d2fb_dphi2(phi::Real, p::Parameters)
    phi_c = clip01(phi, p.eps)
    (1.0 / p.nb) * (1.0 / phi_c + 1.0 / (1 - phi_c)) - 2.0 * p.chi_b
end

# Scalar psi_linear: returns (psi, dpsi, lam).
function psi_linear(z::Real, p::Parameters)
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(alpha / p.kappa)
    psi  = -h0 * p.R * exp(-lam * z) / ((z + p.R) * p.kappa * (lam + 1 / p.R))
    dpsi = psi * (-lam - 1 / (z + p.R))
    return psi, dpsi, lam
end

# Broadcast-aware: applied to an array, returns (psi_array, dpsi_array, lam).
function psi_linear(z::AbstractArray, p::Parameters)
    alpha = d2fb_dphi2(p.phi_inf, p)
    h0  = p.omega_b
    lam = sqrt(alpha / p.kappa)
    psi  = @. -h0 * p.R * exp(-lam * z) / ((z + p.R) * p.kappa * (lam + 1 / p.R))
    dpsi = @. psi * (-lam - 1 / (z + p.R))
    return psi, dpsi, lam
end

compute_mu_inf(p::Parameters) = p.nb * dfb_dphi(p.phi_inf, p)

function compute_spinodal(p::Parameters)
    p.chi_b <= 0 && return (nothing, nothing)
    discriminant = 1 - 2 / (p.chi_b * p.nb)
    discriminant < 0 && return (nothing, nothing)
    sqrt_d = sqrt(discriminant)
    return (0.5 * (1 - sqrt_d), 0.5 * (1 + sqrt_d))
end

# ------------------------------------------------------------------
# 1D isolated-sphere BVP (for the single-sphere-limit comparison test)
# ------------------------------------------------------------------

struct IsolatedSphere1D
    z::Vector{Float64}
    phi::Vector{Float64}
    dphi::Vector{Float64}
    p::Parameters
    success::Bool
end

@inline function _linterp(zg::AbstractVector, yg::AbstractVector,
                          zq::Real, dz::Real, z0::Real)
    n = length(zg)
    t = (zq - z0) / dz
    i = clamp(Int(floor(t)), 0, n - 2)
    a = t - i
    yg[i + 1] * (1 - a) + yg[i + 2] * a
end

function (s::IsolatedSphere1D)(zq::Real)
    z = s.z
    dz = z[2] - z[1]
    zc = clamp(zq, z[1], z[end])
    _linterp(z, s.phi, zc, dz, z[1])
end

function (s::IsolatedSphere1D)(zq::AbstractVector)
    [s(z) for z in zq]
end

# Robin surface slope phi'(0).  The plain BC is -kappa*phi'(0) + omega_b = 0;
# with bc_reg > 0 a regularization term is added (mirrors make_bc in
# viz_util/bvp_solver.py, line 360):
#     -kappa*phi'(0) + omega_b + bc_reg * log(phi_s / (1 - phi_s)) = 0.
# As phi_s -> 1 the log term diverges and opposes the attractive omega_b,
# pinning a unique *finite* adsorption branch near a surface phase transition
# (where the bare BC is bistable: a thin film and a runaway phi_s -> 1).
@inline function _surf_slope(phi_s::Real, omega_b::Real, p::Parameters, bc_reg::Real)
    bc_reg == 0.0 && return omega_b / p.kappa
    ps = clip01(phi_s, p.eps)
    return (omega_b + bc_reg * log(ps / (1.0 - ps))) / p.kappa
end

# d/dphi_s of _surf_slope (for the analytic Jacobian).
@inline function _surf_slope_deriv(phi_s::Real, p::Parameters, bc_reg::Real)
    bc_reg == 0.0 && return 0.0
    ps = clip01(phi_s, p.eps)
    (p.eps < phi_s < 1.0 - p.eps) || return 0.0   # in the clamp the slope is flat
    return bc_reg / (p.kappa * ps * (1.0 - ps))
end

function solve_isolated_sphere(p::Parameters;
                               L_iso::Real = 10.0,
                               N::Int      = 400,
                               f_tol::Real = 1e-8,
                               maxiter::Int = 200,
                               bc_reg::Real = 0.0,
                               omega_continuation::Bool = true)
    @assert isfinite(p.R) "solve_isolated_sphere requires finite p.R"
    z  = collect(range(0.0, L_iso; length = N))
    dz = z[2] - z[1]
    mu_inf  = compute_mu_inf(p)
    Nu = N - 1                                # unknowns at z[1:N-1]
    phi_far = p.phi_inf                       # Dirichlet at z[N]

    # Residual and analytic tridiagonal Jacobian at effective surface field ob.
    function residual!(F, u, ob)
        @inbounds for i in 1:Nu
            phi_i = u[i]
            if i == 1
                # Robin BC ghost: phi(-dz) = phi(+dz) - 2*dz*sk.
                sk    = _surf_slope(phi_i, ob, p, bc_reg)
                ddphi = (2 * u[2] - 2 * phi_i - 2 * dz * sk) / dz^2
                lap   = ddphi + (2.0 / p.R) * sk
            else
                phi_l = u[i - 1]
                phi_r = i == Nu ? phi_far : u[i + 1]
                ddphi = (phi_r - 2 * phi_i + phi_l) / dz^2
                dphi  = (phi_r - phi_l) / (2 * dz)
                rr    = p.R + z[i]
                lap   = ddphi + (2.0 / rr) * dphi
            end
            F[i] = p.kappa * lap - (dfb_dphi(phi_i, p) - mu_inf / p.nb)
        end
        return nothing
    end
    function jacobian!(dl, dd, du, u, ob)
        @inbounds for i in 1:Nu
            if i == 1
                sk1   = _surf_slope_deriv(u[1], p, bc_reg)
                dd[i] = p.kappa * (-2.0 / dz^2) - d2fb_dphi2(u[i], p) +
                        p.kappa * sk1 * (-2.0 / dz + 2.0 / p.R)
                du[i] = p.kappa * (2.0 / dz^2)
                dl[i] = 0.0
            else
                rr    = p.R + z[i]
                dd[i] = p.kappa * (-2.0 / dz^2) - d2fb_dphi2(u[i], p)
                dl[i] = p.kappa * (1.0 / dz^2 - (1.0 / rr) / dz)
                du[i] = i == Nu ? 0.0 : p.kappa * (1.0 / dz^2 + (1.0 / rr) / dz)
            end
        end
        return nothing
    end

    # Damped (backtracking line search) Newton with the analytic tridiagonal
    # Jacobian.  Returns (u, |F|).  Robust on the stiff prewetting branch where
    # the bare dense Newton stalls.
    function newton(u0, ob)
        u  = copy(u0)
        F  = Vector{Float64}(undef, Nu)
        dl = Vector{Float64}(undef, Nu)
        dd = Vector{Float64}(undef, Nu)
        du = Vector{Float64}(undef, Nu)
        ut = Vector{Float64}(undef, Nu)
        residual!(F, u, ob); nrm = maximum(abs, F)
        for _ in 1:maxiter
            nrm < f_tol && break
            jacobian!(dl, dd, du, u, ob)
            δ = Tridiagonal(dl[2:Nu], dd[1:Nu], du[1:Nu - 1]) \ F
            λ = 1.0; accepted = false
            for _ in 1:40
                @inbounds @. ut = clamp(u - λ * δ, p.eps, 1.0 - p.eps)
                residual!(F, ut, ob); newnrm = maximum(abs, F)
                if newnrm < nrm
                    @inbounds @. u = ut; nrm = newnrm; accepted = true; break
                end
                λ *= 0.5
            end
            accepted || break
        end
        return u, nrm
    end

    # Initial guess: exponential decay from a perturbed phi_inf.
    d2  = max(d2fb_dphi2(p.phi_inf, p), 0.05)
    ell = sqrt(p.kappa / d2)
    u0  = Vector{Float64}(undef, Nu)
    @inbounds for i in 1:Nu
        u0[i] = clip01(p.phi_inf - p.omega_b * ell / p.kappa * exp(-z[i] / ell), p.eps)
    end

    u, nrm = newton(u0, p.omega_b)
    # omega_b continuation from the uniform (omega_b = 0) state, tracing the
    # branch continuously when the cold start fails (same strategy as
    # solve_two_spheres, and the branch the 2D solver follows).
    if nrm >= f_tol && omega_continuation && p.omega_b != 0.0
        u_cur = fill(clip01(p.phi_inf, p.eps), Nu)
        n_steps = 8
        for k in 1:n_steps
            u_cur, nrm = newton(u_cur, (k / n_steps) * p.omega_b)
            nrm >= f_tol && break
        end
        nrm < f_tol && (u = u_cur)
    end
    success = nrm < f_tol

    phi      = vcat(u, phi_far)
    dphi     = similar(phi)
    dphi[1]  = _surf_slope(phi[1], p.omega_b, p, bc_reg)
    @inbounds for i in 2:N - 1
        dphi[i] = (phi[i + 1] - phi[i - 1]) / (2 * dz)
    end
    dphi[N] = (phi[N] - phi[N - 1]) / dz

    return IsolatedSphere1D(z, phi, dphi, p, success)
end
