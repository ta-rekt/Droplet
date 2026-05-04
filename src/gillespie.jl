
export buildOffsets!, buildMoves!, computeRate!, deltaH_kmc, decodeSlot, isOccupied,
       fenwick_build!, fenwick_update!, fenwick_prefix, fenwick_total, fenwick_find,
       sampleMove, applyMove!, refreshAffected!, rebuildQ!,
       kmc_step!, kmc_run!

# direction encoding for KMC hops (axis ∈ 1..3, sign ∈ {+1,-1})
const KMC_DIRS = ((Int32(1), Int32(+1)), (Int32(1), Int32(-1)),
                  (Int32(2), Int32(+1)), (Int32(2), Int32(-1)),
                  (Int32(3), Int32(+1)), (Int32(3), Int32(-1)))


############# Phase 1: rate primitives (BulkMonomer-only) #############
# MARK: -Phase1

# fills species_offset and nslots (6 slots per particle, per species)
function buildOffsets!(sys::KMCSystem)
    n = length(sys.slist)
    resize!(sys.species_offset, n)
    cum = 0
    for (i, s) in enumerate(sys.slist)
        sys.species_offset[i] = cum
        cum += 6 * Int(s.length)
    end
    sys.nslots = cum
    return sys
end

# slot id → (speciesIdx, particleId, dir)
@inline function decodeSlot(sys::KMCSystem, slot::Int)
    sIdx = 1
    @inbounds for i in length(sys.species_offset):-1:1
        if slot > sys.species_offset[i]
            sIdx = i
            break
        end
    end
    rel = slot - sys.species_offset[sIdx] - 1
    return sIdx, div(rel, 6) + 1, mod(rel, 6) + 1
end

# pure ΔH for a unit hop of BulkMonomer b from src to dst (no struct mutation
# beyond nnp's transient pos arithmetic, which restores src/dst on exit)
function deltaH_kmc(b::BulkMonomer, src::Vector{Int32}, dst::Vector{Int32})
    @assert b.couplingMat === nothing "deltaH_kmc: couplings not supported in v1"
    dh = 0.0
    for (sp, idx) in b.intSpecies
        dn = nnp(b, sp, dst) - nnp(b, sp, src)
        # the moving particle still sits at src, so it is counted as a self
        # neighbor of dst (src ∈ NN(dst) for any unit hop). Subtract 1.
        if b.intIdx == idx
            dn -= 1
        end
        dh += b.intMat[b.intIdx, idx] * dn
    end
    return dh
end

# is `pos` occupied by any species in slist?
@inline function isOccupied(sys::KMCSystem, pos::Vector{Int32})
    @inbounds for s in sys.slist
        if s.latS[pos[1], pos[2], pos[3]] > 0
            return true
        end
    end
    return false
end

# compute rate for a single slot (no side effects on sys.rates / sys.Q)
function computeRate!(sys::KMCSystem, slot::Int)
    sIdx, particleId, dir = decodeSlot(sys, slot)
    s = sys.slist[sIdx]
    @assert s isa BulkMonomer "Phase 1 supports BulkMonomer only"

    pos_arr = s.pos[Int32(particleId)]
    sys.src[1] = pos_arr[1]; sys.src[2] = pos_arr[2]; sys.src[3] = pos_arr[3]

    axis, sgn = KMC_DIRS[dir]
    sys.dst[1] = sys.src[1]; sys.dst[2] = sys.src[2]; sys.dst[3] = sys.src[3]
    sys.dst[axis] += sgn
    boxmod!(sys.dst, sys.dst, Int32(1), s.boxlength)

    if isOccupied(sys, sys.dst)
        return 0.0
    end

    dh = deltaH_kmc(s, sys.src, sys.dst)
    return s.k0 * min(1.0, exp(-sys.beta * dh))
end

# recompute every rate from scratch and rebuild the Fenwick tree.
# Used by buildMoves! at startup and as a periodic drift-correction pass.
function rebuildQ!(sys::KMCSystem)
    @assert sys.nslots > 0 "rebuildQ!: call buildMoves! first"
    for slot in 1:sys.nslots
        sys.rates[slot] = computeRate!(sys, slot)
    end
    fenwick_build!(sys.fenwick, sys.rates)
    sys.Q = fenwick_total(sys.fenwick)
    return sys
end

# allocate slot arrays from slist sizes, then fill rates + Fenwick.
function buildMoves!(sys::KMCSystem)
    buildOffsets!(sys)
    resize!(sys.rates, sys.nslots)
    fill!(sys.rates, 0.0)
    rebuildQ!(sys)
    return sys
end


############# Phase 2: Fenwick tree + sampling #############
# MARK: -Phase2

# Build Fenwick (BIT) tree in-place from a vector of values. 1-indexed.
function fenwick_build!(tree::Vector{Float64}, vals::Vector{Float64})
    n = length(vals)
    resize!(tree, n)
    copyto!(tree, vals)
    @inbounds for i in 1:n
        j = i + (i & -i)
        if j <= n
            tree[j] += tree[i]
        end
    end
    return tree
end

# Add `delta` to slot i.
@inline function fenwick_update!(tree::Vector{Float64}, i::Int, delta::Float64)
    n = length(tree)
    @inbounds while i <= n
        tree[i] += delta
        i += i & -i
    end
    return tree
end

# Sum of values 1..i.
@inline function fenwick_prefix(tree::Vector{Float64}, i::Int)
    s = 0.0
    @inbounds while i > 0
        s += tree[i]
        i -= i & -i
    end
    return s
end

@inline fenwick_total(tree::Vector{Float64}) = fenwick_prefix(tree, length(tree))

# Smallest i such that prefix(i) > u. With u ∈ [0, total), this skips zero-rate
# slots: any slot with rate 0 leaves the prefix unchanged, so prefix(i) > u
# selects only slots that contribute strictly positive weight.
function fenwick_find(tree::Vector{Float64}, u::Float64)
    n = length(tree)
    bit = 1
    while bit <= n
        bit <<= 1
    end
    bit >>= 1
    i = 0
    @inbounds while bit > 0
        nxt = i + bit
        if nxt <= n && tree[nxt] <= u
            i = nxt
            u -= tree[nxt]
        end
        bit >>= 1
    end
    return i + 1
end

# Sample a slot proportional to its rate. Caller ensures sys.Q > 0.
@inline function sampleMove(sys::KMCSystem)
    @assert sys.Q > 0.0 "sampleMove: total rate is zero"
    return fenwick_find(sys.fenwick, rand() * sys.Q)
end


############# Phase 3: applyMove! + incremental updates #############
# MARK: -Phase3

# inverse direction (KMC_DIRS[d] points the opposite way to KMC_DIRS[reverse_dir(d)])
@inline reverse_dir(d::Integer) = ((d - 1) ⊻ 1) + 1

# refresh every slot whose source or destination is in NN1(p_old) ∪ NN1(p_new),
# then reconcile sys.Q and the Fenwick tree by the per-slot deltas.
# IMPORTANT: p_old/p_new must not alias sys.src or sys.dst once the refresh
# loop starts, since computeRate! overwrites those buffers. The collection
# phase (which reads p_old/p_new) finishes before any computeRate! call.
function refreshAffected!(sys::KMCSystem, p_old::Vector{Int32}, p_new::Vector{Int32})
    empty!(sys.refresh_set)
    _collectAround!(sys, p_old)
    _collectAround!(sys, p_new)

    for slot in sys.refresh_set
        old_rate = sys.rates[slot]
        new_rate = computeRate!(sys, slot)
        if new_rate != old_rate
            sys.rates[slot] = new_rate
            delta = new_rate - old_rate
            fenwick_update!(sys.fenwick, slot, delta)
            sys.Q += delta
        end
    end
    return sys
end

# enumerate NN1(center) (center plus its 6 face neighbors) and gather the slots
# whose source or destination lands at any of those sites.
function _collectAround!(sys::KMCSystem, center::Vector{Int32})
    boxlength = first(sys.slist).boxlength
    _addSiteAffected!(sys, center)

    nb = sys.nbuf
    for d in 1:6
        axis, sgn = KMC_DIRS[d]
        nb[1] = center[1]; nb[2] = center[2]; nb[3] = center[3]
        nb[axis] += sgn
        boxmod!(nb, nb, Int32(1), boxlength)
        _addSiteAffected!(sys, nb)
    end
    return sys
end

# at site x: any particle at x has its 6 source-moves added; any particle at a
# face-neighbor of x has its move-toward-x added (covers dest-at-x cases).
function _addSiteAffected!(sys::KMCSystem, x::Vector{Int32})
    boxlength = first(sys.slist).boxlength

    # source-at-x
    for (sIdx, s) in enumerate(sys.slist)
        pid = s.latS[x[1], x[2], x[3]]
        if pid > 0
            base = sys.species_offset[sIdx] + (Int(pid) - 1) * 6
            for d in 1:6
                push!(sys.refresh_set, base + d)
            end
        end
    end

    # dest-at-x: walk NN(x); particles there have a move toward x with dir = reverse_dir(d)
    yp = sys.ybuf
    for d in 1:6
        axis, sgn = KMC_DIRS[d]
        yp[1] = x[1]; yp[2] = x[2]; yp[3] = x[3]
        yp[axis] += sgn
        boxmod!(yp, yp, Int32(1), boxlength)
        rev = reverse_dir(d)
        for (sIdx, s) in enumerate(sys.slist)
            pid = s.latS[yp[1], yp[2], yp[3]]
            if pid > 0
                push!(sys.refresh_set,
                      sys.species_offset[sIdx] + (Int(pid) - 1) * 6 + rev)
            end
        end
    end
    return sys
end

# Apply the hop encoded by `slot`: advance simulation time using the OLD Q,
# update lattice + pos, then refresh affected rates incrementally.
function applyMove!(sys::KMCSystem, slot::Int)
    sIdx, particleId, dir = decodeSlot(sys, slot)
    s = sys.slist[sIdx]
    @assert s isa BulkMonomer "Phase 3 supports BulkMonomer only"

    # source and destination — held in sys.src/sys.dst, then read by refreshAffected!
    # during its collection pass (before any computeRate! call clobbers them).
    pos_arr = s.pos[Int32(particleId)]
    sys.src[1] = pos_arr[1]; sys.src[2] = pos_arr[2]; sys.src[3] = pos_arr[3]

    axis, sgn = KMC_DIRS[dir]
    sys.dst[1] = sys.src[1]; sys.dst[2] = sys.src[2]; sys.dst[3] = sys.src[3]
    sys.dst[axis] += sgn
    boxmod!(sys.dst, sys.dst, Int32(1), s.boxlength)

    # Gillespie time advance using the rate distribution that produced this draw
    Q_old = sys.Q
    @assert Q_old > 0.0 "applyMove!: Q == 0"
    sys.t_sim += -log(rand()) / Q_old
    sys.n_events += 1

    # commit the move on lattice + pos
    s.latS[sys.src[1], sys.src[2], sys.src[3]] = Int32(0)
    s.latS[sys.dst[1], sys.dst[2], sys.dst[3]] = Int32(particleId)
    pos_arr[1] = sys.dst[1]; pos_arr[2] = sys.dst[2]; pos_arr[3] = sys.dst[3]

    refreshAffected!(sys, sys.src, sys.dst)
    return sys
end


############# Phase 4: main loop + sim-time sampling #############
# MARK: -Phase4

# one Gillespie event: sample → applyMove! (which advances t_sim and refreshes)
@inline function kmc_step!(sys::KMCSystem)
    applyMove!(sys, sampleMove(sys))
    return sys
end

# Run KMC until any of the caps fire. Sampling is on a fixed grid of simulation
# time t_burn_in, t_burn_in + dt_sample, t_burn_in + 2·dt_sample, ...; whenever
# sys.t_sim crosses the next grid point, sample_callback(sys) is invoked once.
# t_max_wall is in minutes (matches target!() convention).
function kmc_run!(sys::KMCSystem;
                  t_max_wall::Float64=Inf,
                  t_max_sim::Float64=Inf,
                  t_burn_in::Float64=0.0,
                  dt_sample::Float64=Inf,
                  sample_callback=nothing,
                  rebuild_every::Int=1_000_000)
    t_start = time_ns()
    t_next_sample = t_burn_in
    next_rebuild = sys.n_events + rebuild_every
    n_samples = 0

    while sys.t_sim < t_max_sim
        if (time_ns() - t_start) / 6e10 >= t_max_wall
            break
        end
        if sys.Q <= 0.0
            break
        end

        kmc_step!(sys)

        if sys.n_events >= next_rebuild
            rebuildQ!(sys)
            next_rebuild = sys.n_events + rebuild_every
        end

        if sample_callback !== nothing && sys.t_sim >= t_next_sample
            sample_callback(sys)
            n_samples += 1
            t_next_sample += dt_sample
        end
    end

    return n_samples
end
