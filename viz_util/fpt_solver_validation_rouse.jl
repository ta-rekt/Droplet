# First-passage-time OUT OF the surface-adsorption well, in the combined
# attractive-adsorption + Rouse-entropy potential.  Three independent estimates
# of the same escape-time distribution are compared:
#
#   1. backward-FPE solver (src/fpt_solver.jl) on the reduced 1D coordinate,
#   2. direct 1D Euler–Maruyama Langevin in the SAME reduced potential U_tot,
#   3. the FULL Rouse polymer (src/rouse.jl): a 151-bead chain whose two ends
#      feel the adsorption force, with escape times read from the A-sojourn
#      (bound-state dwell) machinery.
#
# The reduced coordinate x is the end-to-end (bead-to-bead) distance.  U_tot
# already contains the Rouse entropy (-2 log x is the radial Jacobian, plus the
# Gaussian spring), so the 1D solver with a flat operator is the exact radial
# Smoluchowski equation for that coordinate.  Its diffusion coefficient is the
# relative-coordinate value D = 2 kBT/gamma (two end beads).  The full Rouse
# chain applies only the adsorption force U_ads between the ends — its spring and
# entropy emerge from the chain, so this tests the 1D reduction against the
# genuine many-bead dynamics.
#
#   julia --project=@. viz_util/fpt_solver_validation_rouse.jl
#
# Writes PNG + PDF to media/figures/<today>/.

ENV["GKSwstype"] = "100"

using Droplet
using Plots
using Random
using Statistics
using Printf
using Dates
import JSON

# --- Combined potential ------------------------------------------------------
p = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.25, kappa = 1.0,
               phi_inf = 0.1, eps = 1e-6, R = 5.0)     # omega_b weakened to ~7 kBT
d2f_phi_inf = d2fb_dphi2(p.phi_inf, p)
h0  = p.omega_b
lam = sqrt(d2f_phi_inf / p.kappa)
kuhnLength = 1.0
beadSep    = 150

U_ads(x) = -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
            (p.kappa * (lam + 1 / p.R) * x)
U_tot(x) = U_ads(x) - 2 * log(x) + 3 * x^2 / (2 * beadSep * kuhnLength^2)

h_fd = 1e-6
U_tot_prime(x) = -(U_tot(x + h_fd) - U_tot(x - h_fd)) / (2 * h_fd)
U_ads_force(x) = -(U_ads(x + h_fd) - U_ads(x - h_fd)) / (2 * h_fd)  # -dU_ads/dr

# --- Physical constants (shared by all three methods) ------------------------
kBT      = 1.0
friction = 1.0
D_bead   = kBT / friction
D_eff    = 2 * D_bead          # relative coordinate of the two end beads

# --- Locate the ridge (absorbing point) --------------------------------------
xs = collect(range(p.R, 25.0; length = 8000))
Us = U_tot.(xs)
ridge_idx = findfirst(i -> Us[i] > Us[i-1] && Us[i] > Us[i+1], 2:length(xs)-1)
x_ridge = ridge_idx === nothing ? xs[argmax(Us)] : xs[ridge_idx + 1]
barrier = U_tot(x_ridge) - U_tot(p.R)
@printf("lam=%.3f  ridge x=%.3f  barrier=%.2f kBT\n", lam, x_ridge, barrier)

# --- 1. Backward-FPE solver (IC = delta at R, reflecting at R, absorbing at ridge)
xg = collect(range(p.R, x_ridge; length = 2000))
Ug = U_tot.(xg)
p0 = zeros(length(xg)); p0[1] = 1.0                    # start at contact x = R
prob = FPTProblem(xg, Ug; D = D_eff, kBT = kBT,
                  left_bc = :reflecting, right_bc = :absorbing, p0 = p0)
res  = solve_fpt(prob)

# --- 2. Direct 1D Langevin in U_tot (IC = R) ---------------------------------
# `force` is -dU/dx (same convention as U_tot_prime and the Rouse pair force);
# the overdamped drift is D*β*force*dt.
function bd_fpt_from_R(force, xmin, xabs; D, kBT, dt, nsamples, seed)
    rng = MersenneTwister(seed)
    β = 1.0 / kBT
    c1 = D * β * dt
    c2 = sqrt(2 * D * dt)
    xa = xabs - 0.5826 * sqrt(2 * D * dt)   # Milstein absorbing-boundary correction
    times = Float64[]
    maxstep = 20_000_000
    for _ in 1:nsamples
        xc = xmin
        for step in 1:maxstep
            xc += c1 * force(xc) + c2 * randn(rng)
            xc <= xmin && (xc = 2 * xmin - xc)
            xc >= xa && (push!(times, step * dt); break)
        end
    end
    return times
end

dt_bd = 2e-4            # small dt: the well is steep, EM biases the MFPT otherwise
t_bd  = bd_fpt_from_R(U_tot_prime, p.R, x_ridge;
                      D = D_eff, kBT = kBT, dt = dt_bd, nsamples = 10_000, seed = 7)

# --- 3. Full Rouse polymer (151 beads); escape times from A-sojourns ---------
a_bound = p.R + 0.4        # bound-state cutoff (r_min = R is the reflecting wall)
b_open  = x_ridge          # ridge = escape threshold
dt_r    = 5e-4
nsteps  = 16_000_000
r = RousePolymer(; boltzmannConst = 1.0, temperature = kBT, friction = friction,
                 kuhnLength = kuhnLength, r_min = p.R,
                 nbeads = Int32(beadSep + 1), nsteps = Int32(nsteps), dt = dt_r,
                 trackedBead1 = 1, trackedBead2 = 1 + beadSep,
                 attractiveForce = U_ads_force, saveMode = :sojourn_only,
                 a = a_bound, b = b_open, seed = 20260701)
init!(r); simulate!(r)
t_rouse = Float64.(r.A_intervals) .* dt_r      # bound-state escape times

@printf("MFPT   solver=%.3f   1D-Langevin=%.3f (n=%d)   full-Rouse=%.3f (n=%d escapes)\n",
        res.mfpt, mean(t_bd), length(t_bd),
        isempty(t_rouse) ? NaN : mean(t_rouse), length(t_rouse))

# --- Figure ------------------------------------------------------------------
# potential panel
xplot = collect(range(p.R, min(2 * x_ridge, 14.0); length = 600))
pp = plot(xplot, U_tot.(xplot); lw = 2.5, color = :black, legend = false,
          xlabel = "bead-to-bead distance x", ylabel = "U_tot(x) / k_BT",
          title = @sprintf("Combined potential (barrier %.1f k_BT)", barrier))
vline!(pp, [p.R];     color = :dodgerblue, ls = :dash, lw = 2)   # reflecting (contact)
vline!(pp, [x_ridge]; color = :crimson,    ls = :solid, lw = 2)  # absorbing (ridge)
annotate!(pp, p.R + 0.15, U_tot(p.R) + 0.5, text("well (R)", 8, :left, :dodgerblue))
annotate!(pp, x_ridge - 0.15, U_tot(x_ridge) - 1.2, text("ridge", 8, :right, :crimson))

# FPT panel: solver vs 1D Langevin vs full Rouse (semilog-y).  The full-Rouse
# sample is smaller, so it is drawn as a step outline rather than filled bars.
thi  = quantile([t_bd; t_rouse], 0.99)
keep = res.t .<= thi
ymax = maximum(res.pdf[keep])
edges = range(0, thi; length = 45)
pf = histogram(t_bd; bins = edges, normalize = :pdf, color = :dodgerblue,
               alpha = 0.35, linecolor = :white, label = "1D Langevin (n=$(length(t_bd)))",
               yscale = :log10, ylims = (ymax / 200, ymax * 1.6), xlims = (0, thi),
               xlabel = "escape time t", ylabel = "PDF  f(t)",
               title = @sprintf("Escape out of the well:  MFPT solver %.2f / BD %.2f / Rouse %.2f",
                                res.mfpt, mean(t_bd), isempty(t_rouse) ? NaN : mean(t_rouse)))
stephist!(pf, t_rouse; bins = edges, normalize = :pdf, color = :seagreen, lw = 2.5,
          label = "full Rouse (n=$(length(t_rouse)))")
plot!(pf, res.t[keep], res.pdf[keep]; lw = 3, color = :crimson, label = "backward-FPE solver")
vline!(pf, [res.mfpt]; color = :crimson, ls = :dot, lw = 2, label = "MFPT (solver)")

fig = plot(pp, pf; layout = grid(1, 2, widths = [0.36, 0.64]),
           size = (1600, 620), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
           dpi = 150,
           plot_title = "FPT out of the adsorption well: 1D solver vs 1D Langevin vs full Rouse chain")

outdir = joinpath("media", "figures", Dates.format(today(), "yyyymmdd"))
mkpath(outdir)
png = joinpath(outdir, "fpt_solver_validation_rouse.png")
pdf = joinpath(outdir, "fpt_solver_validation_rouse.pdf")
savefig(fig, png); savefig(fig, pdf)
println("saved -> ", png)
println("saved -> ", pdf)

# dump the samples so the figure can be re-plotted without re-running the sim
open(joinpath(outdir, "fpt_solver_validation_rouse.json"), "w") do io
    JSON.print(io, Dict(
        "omega_b" => p.omega_b, "barrier" => barrier, "x_ridge" => x_ridge,
        "D_eff" => D_eff, "beadSep" => beadSep, "dt_bd" => dt_bd, "dt_rouse" => dt_r,
        "mfpt_solver" => res.mfpt, "mfpt_bd" => mean(t_bd),
        "mfpt_rouse" => isempty(t_rouse) ? NaN : mean(t_rouse),
        "solver_t" => res.t, "solver_pdf" => res.pdf,
        "t_bd" => t_bd, "t_rouse" => t_rouse))
end