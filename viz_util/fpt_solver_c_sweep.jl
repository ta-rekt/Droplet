# Concentration sweep of the first-passage-time distributions in the combined
# attractive-adsorption + Rouse-entropy potential (the reduced 1D end-to-end
# coordinate x, cf. fpt_solver_validation_rouse.jl).
#
# The bulk volume fraction phi_inf is the only swept parameter: it sets the
# screening length lam = sqrt(d2fb_dphi2(phi_inf) / kappa), which controls the
# range/depth of the adsorption well U_ads.  Higher phi_inf => smaller lam =>
# deeper, longer-ranged well => higher escape barrier.
#
# For every phi_inf the backward-FPE solver (src/fpt_solver.jl) is run TWICE on
# the domain [R, x_ridge]:
#
#   1. escape  ("out of the well"): start at contact x = R (reflecting wall),
#      absorbing at the ridge  -> time to climb OUT of the well and reach x_ridge.
#   2. return  ("into the well"):   start at the ridge on the right (reflecting
#      wall), absorbing at x = R  -> time to reach x = R starting from the right.
#
# Both densities, their means (MFPTs), the potential on the grid, and every
# fixed/varied parameter are dumped to data/<today>/fpt_solver_c_sweep/.
#
#   julia --project=@. viz_util/fpt_solver_c_sweep.jl

ENV["GKSwstype"] = "100"

using Droplet
using Plots
using Statistics
using Printf
using Dates
import JSON

# --- Fixed parameters (only phi_inf is varied) -------------------------------
phi_infs   = collect(range(0.001, stop = 0.05, length = 50))
p_base     = Parameters(nb = 1.0, chi_b = 2.5, omega_b = -0.25, kappa = 1.0,
                        phi_inf = 0.1, eps = 1e-6, R = 5.0)
h0         = p_base.omega_b        # adsorption strength, fixed across the sweep
kuhnLength = 1.0
beadSep    = 150
kBT        = 1.0
friction   = 1.0
D_bead     = kBT / friction
D_eff      = 2 * D_bead            # relative coord of the two end beads
Ngrid      = 1500                  # grid points on [R, x_ridge]

# combined potential for a given screening length lam (h0, beadSep etc. fixed)
make_U_tot(p, lam) =
    x -> -(4π * p.R^3 * h0^2 * exp(-lam * (x - p.R))) /
          (p.kappa * (lam + 1 / p.R) * x) -              # U_ads
         2 * log(x) +                                     # radial Jacobian
         3 * x^2 / (2 * beadSep * kuhnLength^2)           # Rouse spring

# --- Sweep -------------------------------------------------------------------
records = Vector{Dict{String,Any}}()
for k in eachindex(phi_infs[1:end-2])
    φ = phi_infs[k]
    φn = phi_infs[k + 1]
    φnn = phi_infs[k + 2]

    # fix t_grid, calculate second derivative of pdf wrsp to φ, integrate and record it

    p   = Droplet.with(p_base; phi_inf = φ)
    lam = sqrt(d2fb_dphi2(φ, p) / p.kappa)
    U_tot = make_U_tot(p, lam)

    # locate the ridge (barrier top = absorbing point for escape)
    xs = collect(range(p.R, 25.0; length = 8000))
    Us = U_tot.(xs)
    ridge_idx = findfirst(i -> Us[i] > Us[i-1] && Us[i] > Us[i+1], 2:length(xs)-1)
    x_ridge = ridge_idx === nothing ? xs[argmax(Us)] : xs[ridge_idx + 1]
    barrier = U_tot(x_ridge) - U_tot(p.R)

    # potential on the [R, ridge] grid shared by both directions
    xg_esc = collect(range(p.R, x_ridge; length = Ngrid))
    Ug_esc = U_tot.(xg_esc)

    # 1. escape: R -> ridge (start at contact, reflecting at R, absorbing at ridge)
    p0_esc = zeros(Ngrid); p0_esc[1] = 1.0
    prob_esc = FPTProblem(xg_esc, Ug_esc; D = D_eff, kBT = kBT, t=t_esc_grid,
                          left_bc = :reflecting, right_bc = :absorbing, p0 = p0_esc)
    res_esc = solve_fpt(prob_esc)

    # 2. return: ridge -> R (start at ridge, absorbing at R, reflecting at infty)
    xg_ret = collect(range(p.R, 4*kuhnLength*sqrt(beadSep); length = Ngrid))
    Ug_ret = U_tot.(xg_ret)
    p0_ret = zeros(Ngrid); p0_ret[ridge_idx] = 1.0
    prob_ret = FPTProblem(xg_ret, Ug_ret; D = D_eff, kBT = kBT, t=t_ret_grid,
                          left_bc = :absorbing, right_bc = :reflecting, p0 = p0_ret)
    res_ret = solve_fpt(prob_ret)

    push!(records, Dict(
        "phi_inf"     => φ,
        "lam"         => lam,
        "x_ridge"     => x_ridge,
        "barrier"     => barrier,
        "well_depth"  => U_tot(p.R),
        "x_grid_esc"  => xg_esc,
        "U_tot_esc"   => Ug_esc,
        "x_grid_ret"  => xg_ret,
        "U_tot_ret"   => Ug_ret,
        "escape_t"    => res_esc.t, "escape_pdf" => res_esc.pdf, "escape_mfpt" => res_esc.mfpt,
        "return_t"    => res_ret.t, "return_pdf" => res_ret.pdf, "return_mfpt" => res_ret.mfpt))

    @printf("phi=%.4f  lam=%5.2f  ridge=%.3f  barrier=%.2f kBT  MFPT_esc=%.4g  MFPT_ret=%.4g\n",
            φ, lam, x_ridge, barrier, res_esc.mfpt, res_ret.mfpt)
end

# --- Dump everything to data/<today>/fpt_solver_c_sweep/ ----------------------
outdir = joinpath("data", Dates.format(today(), "yyyymmdd"), "fpt_solver_c_sweep")
mkpath(outdir)
jsonpath = joinpath(outdir, "fpt_solver_c_sweep.json")
open(jsonpath, "w") do io
    JSON.print(io, Dict(
        "fixed" => Dict(
            "nb" => p_base.nb, "chi_b" => p_base.chi_b, "omega_b" => p_base.omega_b,
            "kappa" => p_base.kappa, "R" => p_base.R, "eps" => p_base.eps,
            "h0" => h0, "kBT" => kBT, "friction" => friction, "D_eff" => D_eff,
            "beadSep" => beadSep, "kuhnLength" => kuhnLength, "Ngrid" => Ngrid,
            "escape_setup" => "IC=delta at R; reflecting at R, absorbing at ridge (time R->ridge)",
            "return_setup" => "IC=delta at ridge; absorbing at R, reflecting at ridge (time ridge->R)"),
        "phi_infs" => phi_infs,
        "per_phi"  => records))
end
println("saved -> ", jsonpath)

# --- Overview figure: MFPTs and barrier vs concentration ---------------------
mfpt_esc = [r["escape_mfpt"] for r in records]
mfpt_ret = [r["return_mfpt"] for r in records]
barriers = [r["barrier"]     for r in records]

pm = plot(phi_infs, mfpt_esc; lw = 2.5, color = :crimson, yscale = :log10,
          label = "escape  (R -> ridge)", marker = :circle, ms = 3,
          xlabel = "bulk volume fraction phi_inf", ylabel = "MFPT",
          title = "First-passage times vs concentration", legend = :topleft)
plot!(pm, phi_infs, mfpt_ret; lw = 2.5, color = :seagreen, marker = :circle, ms = 3,
      label = "return  (ridge -> R)")
pb = plot(phi_infs, barriers; lw = 2.5, color = :black, marker = :circle, ms = 3,
          xlabel = "bulk volume fraction phi_inf", ylabel = "barrier / k_BT",
          title = "Escape barrier vs concentration", legend = false)
fig = plot(pm, pb; layout = (1, 2), size = (1400, 550), dpi = 150,
           left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
pngpath = joinpath(outdir, "fpt_solver_c_sweep.png")
savefig(fig, pngpath)
println("saved -> ", pngpath)
