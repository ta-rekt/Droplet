# Flory-Huggins bulk thermodynamics + 1D isolated-sphere BVP solver.
#
# Mirrors viz_util/bvp_solver.py: provides the Parameters dataclass and the
# scalar thermodynamic functions (fb, dfb_dphi, d2fb_dphi2, psi_linear,
# compute_mu_inf, compute_spinodal) that the 2D two-sphere code calls into.
#
# Only the minimum 1D solver needed by the two-sphere test suite is ported
# (the test's "single-sphere limit" comparison).  The full robust/shooting/
# continuation pipeline in the Python file is not reproduced.

Base.@kwdef struct Parameters
    nb::Float64              = 1.0
    chi_b::Float64           = 0.0
    omega_b::Float64         = 0.0
    kappa::Float64           = 1.0
    phi_inf::Float64         = 0.06
    L::Float64               = 10.0
    R::Float64               = 3.0
    eps::Float64             = 1e-10
    conserve_volume::Bool    = false
    parallel_plates::Bool    = false
end

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

function solve_isolated_sphere(p::Parameters;
                               L_iso::Real = 10.0,
                               N::Int      = 400,
                               f_tol::Real = 1e-8,
                               maxiter::Int = 200)
    @assert isfinite(p.R) "solve_isolated_sphere requires finite p.R"
    z  = collect(range(0.0, L_iso; length = N))
    dz = z[2] - z[1]
    mu_inf      = compute_mu_inf(p)
    omega_kappa = p.omega_b / p.kappa
    Nu = N - 1                                # unknowns at z[1:N-1]
    phi_far = p.phi_inf                       # Dirichlet at z[N]

    function residual!(F, u, _)
        @inbounds for i in 1:Nu
            phi_i = u[i]
            if i == 1
                # Robin BC ghost: phi(-dz) = phi(+dz) - 2*dz*omega_kappa,
                # giving ddphi = (2*u[2] - 2*u[1] - 2*dz*omega_kappa)/dz^2.
                ddphi = (2 * u[2] - 2 * phi_i - 2 * dz * omega_kappa) / dz^2
                lap   = ddphi + (2.0 / p.R) * omega_kappa
            else
                phi_l = u[i - 1]
                phi_r = i == Nu ? phi_far : u[i + 1]
                ddphi = (phi_r - 2 * phi_i + phi_l) / dz^2
                dphi  = (phi_r - phi_l) / (2 * dz)
                rr    = p.R + z[i]
                lap   = ddphi + (2.0 / rr) * dphi
            end
            driving = dfb_dphi(phi_i, p) - mu_inf / p.nb
            F[i]    = p.kappa * lap - driving
        end
        return nothing
    end

    # Initial guess: exponential decay from a perturbed phi_inf.
    d2  = max(d2fb_dphi2(p.phi_inf, p), 0.05)
    ell = sqrt(p.kappa / d2)
    u0  = Vector{Float64}(undef, Nu)
    @inbounds for i in 1:Nu
        u0[i] = clip01(p.phi_inf - p.omega_b * ell / p.kappa * exp(-z[i] / ell), p.eps)
    end

    prob = NonlinearProblem(residual!, u0)
    sol  = solve(prob, NewtonRaphson(autodiff = AutoFiniteDiff());
                 abstol = f_tol, maxiters = maxiter)
    success = SciMLBase.successful_retcode(sol)

    phi      = vcat(sol.u, phi_far)
    dphi     = similar(phi)
    dphi[1]  = omega_kappa
    @inbounds for i in 2:N - 1
        dphi[i] = (phi[i + 1] - phi[i - 1]) / (2 * dz)
    end
    dphi[N] = (phi[N] - phi[N - 1]) / dz

    return IsolatedSphere1D(z, phi, dphi, p, success)
end
