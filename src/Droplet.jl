module Droplet

    using Statistics
    using StatsBase
    using Distributions
    using Serialization
    using LinearAlgebra
    using JSON
    using DataStructures
    using OrderedCollections
    using TimerOutputs
    using Printf
    using NonlinearSolve
    using NonlinearSolve.ADTypes: AutoFiniteDiff
    using NonlinearSolve.LinearSolve: KrylovJL_GMRES
    using NonlinearSolve.SciMLBase

    include("types.jl")
    include("interactions.jl")
    include("hamiltonian.jl")
    include("metropolis.jl")
    include("move.jl")
    include("port.jl")
    include("init.jl")
    include("utils.jl")
    include("gillespie.jl")
    include("stats.jl")
    include("ctmc.jl")
    include("generate.jl")
    include("get_residence_time.jl")

    # 2D axisymmetric Flory-Huggins BVP (two spheres) + 1D isolated-sphere
    # reference solver.  Independent of the Metropolis core.
    include("bvp_solver.jl")
    include("bvp_two_spheres.jl")

    export Parameters, with, fb, dfb_dphi, d2fb_dphi2, psi_linear,
           compute_mu_inf, compute_spinodal, is_spherical,
           IsolatedSphere1D, solve_isolated_sphere,
           Geometry, build_geometry, decay_length,
           scatter, scatter!, laplacian, laplacian!,
           make_residual, make_residual_inplace,
           initial_guess, solve_two_spheres,
           interp_phi, cell_fluid_fraction,
           phi_on_sphere, excess_free_energy_two_spheres,
           quadratic_ghost_coefficients, find_clean_stencil

end # module Droplet
