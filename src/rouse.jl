
export init!, simulate!

# Rouse connectivity matrix A: tridiagonal with diag (1,2,2,...,2,1), off-diag -1.
function _build_adjacency_matrix(N::Integer)
    A = zeros(Float64, N, N)
    A[1, 1] = 1.0
    A[N, N] = 1.0
    for i in 2:(N - 1)
        A[i, i] = 2.0
    end
    for i in 1:(N - 1)
        A[i, i + 1] = -1.0
        A[i + 1, i] = -1.0
    end
    return A
end

# Characteristic scales: timeScale = γ/κ, posScale = √(2 kBT/κ); κ = 3 kBT/b².
function _set_scaling_params!(r::RousePolymer)
    r.stiffness = 3 * r.boltzmannConst * r.temperature / r.kuhnLength^2
    r.timeScale = r.friction / r.stiffness
    r.posScale  = sqrt(2 * r.boltzmannConst * r.temperature / r.stiffness)
    r.dtDimless = r.dt / r.timeScale
    return r
end

# Allocate working buffers and trajectory storage.
function _set_simulation_vars!(r::RousePolymer)
    N = Int(r.nbeads)
    r.adjacency = _build_adjacency_matrix(N)
    r.pos       = zeros(Float64, N, 3)
    r.drift     = zeros(Float64, N, 3)
    r.noise     = zeros(Float64, N, 3)
    r.distanceTraj = zeros(Float64, r.nsteps)
    if r.saveMode === :full
        r.posTraj = zeros(Float64, r.nsteps, N, 3)
    else
        r.posTraj = nothing
    end
    return r
end

# Overwrite the N×3 noise buffer in-place from rng.
function _set_noise!(r::RousePolymer, rng::AbstractRNG)
    randn!(rng, r.noise)
    return r.noise
end

# Sample an equilibrium Gaussian chain: bond vectors ~ N(0, b²/3) per axis.
function _sample_equilibrium_chain(rng::AbstractRNG, N::Integer, b::Float64)
    pos = zeros(Float64, N, 3)
    sigma = b / sqrt(3)
    @inbounds for i in 2:N
        pos[i, 1] = pos[i - 1, 1] + sigma * randn(rng)
        pos[i, 2] = pos[i - 1, 2] + sigma * randn(rng)
        pos[i, 3] = pos[i - 1, 3] + sigma * randn(rng)
    end
    return pos
end

# Initialize scales, buffers, initial positions, and the first sample of distanceTraj.
function init!(r::RousePolymer)
    _set_scaling_params!(r)
    _set_simulation_vars!(r)

    rng = r.seed === nothing ? Random.default_rng() : Random.MersenneTwister(r.seed)

    if r.initPos === nothing
        r.initPos = _sample_equilibrium_chain(rng, Int(r.nbeads), r.kuhnLength)
    end
    r.pos .= r.initPos ./ r.posScale

    b1, b2 = Int(r.trackedBead1), Int(r.trackedBead2)
    dx = (r.pos[b1, 1] - r.pos[b2, 1]) * r.posScale
    dy = (r.pos[b1, 2] - r.pos[b2, 2]) * r.posScale
    dz = (r.pos[b1, 3] - r.pos[b2, 3]) * r.posScale
    r.distanceTraj[1] = sqrt(dx * dx + dy * dy + dz * dz)

    if r.saveMode === :full
        @views r.posTraj[1, :, :] .= r.pos .* r.posScale
    end
    return r
end

# Dimensionless Euler-Maruyama:
#   x ← x − (A·x)·dτ + ξ·√dτ,   ξ ~ N(0,1) per bead per dim.
# Distance trajectory is recorded each step; positions only when saveMode == :full.
function simulate!(r::RousePolymer)
    rng = r.seed === nothing ? Random.default_rng() : Random.MersenneTwister(r.seed)
    # consume the same draws that init! used so the seed produces reproducible streams
    # (init! constructs its own rng from the same seed -- safe to re-seed here)

    sdt = sqrt(r.dtDimless)
    ps  = r.posScale
    b1  = Int(r.trackedBead1)
    b2  = Int(r.trackedBead2)
    fullMode = r.saveMode === :full
    f   = r.attractiveForce
    invKappaPos = 1 / (r.stiffness * ps)

    @inbounds for t in 2:r.nsteps
        mul!(r.drift, r.adjacency, r.pos)

        # Add the radial pair force between the two tracked beads to drift.
        # Convention: f(d_phys) = −dU/dr; the dimensionless force on i1 is
        #   F/(κ·posScale) = f(d_phys) · (x_i1 − x_i2) / (κ·posScale·|Δx|).
        # We *subtract* it from drift so the existing update
        #   pos ← pos − drift·dτ + noise·√dτ
        # adds the force contribution with the correct sign.
        if f !== nothing
            dxv  = r.pos[b1, 1] - r.pos[b2, 1]
            dyv  = r.pos[b1, 2] - r.pos[b2, 2]
            dzv  = r.pos[b1, 3] - r.pos[b2, 3]
            dmag = sqrt(dxv*dxv + dyv*dyv + dzv*dzv)
            if dmag > 0
                coeff = -invKappaPos * f(dmag * ps) / dmag
                r.drift[b1, 1] += coeff * dxv
                r.drift[b1, 2] += coeff * dyv
                r.drift[b1, 3] += coeff * dzv
                r.drift[b2, 1] -= coeff * dxv
                r.drift[b2, 2] -= coeff * dyv
                r.drift[b2, 3] -= coeff * dzv
            end
        end

        randn!(rng, r.noise)
        @. r.pos = r.pos - r.drift * r.dtDimless + r.noise * sdt

        dx = (r.pos[b1, 1] - r.pos[b2, 1]) * ps
        dy = (r.pos[b1, 2] - r.pos[b2, 2]) * ps
        dz = (r.pos[b1, 3] - r.pos[b2, 3]) * ps
        r.distanceTraj[t] = sqrt(dx * dx + dy * dy + dz * dz)

        if fullMode
            @views r.posTraj[t, :, :] .= r.pos .* ps
        end
    end
    return r
end
