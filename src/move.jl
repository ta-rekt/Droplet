# ensures positions are periodic within the lattice's bounding box
export boxmod!
function boxmod!(v::Union{SubArray{Int32}, Vector{Int32}}, 
                 w::Union{SubArray{Int32}, Vector{Int32}}, 
                 lb::Int32, ub::Int32)

    for i in eachindex(v)
        if v[i] > ub
            w[i] = v[i] - (ub - lb + 1)
        elseif v[i] < lb
            w[i] = v[i] + (ub - lb + 1)
        else
            w[i] = v[i]
        end
    end
    
    return w
end

@inline function coords_equal(a, b)
    return a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
end

export overlap
function overlap(s::Species)
    return ((s.latS[s.newPosi[1], s.newPosi[2], s.newPosi[3]] == 1 && 
             !coords_equal(s.bprev, s.newPosi)) ||
            (s.latS[s.newPosj[1], s.newPosj[2], s.newPosj[3]] == 1 && 
             !coords_equal(s.bnext, s.newPosj)))
end

export addBond!
function addBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi != s.bprev) && (s.newPosj != s.bnext) begin
        "Move is not a bond addition"
    end

    # set new lattice sites to 1
    s.latS[s.newPosi...] = 1
    s.latS[s.newPosj...] = 1

    # Find next available keys
    next_key1 = KEYCOUNT[] + 1
    next_key2 = KEYCOUNT[] + 2
    while haskey(s.pos, next_key1) || haskey(s.pos, next_key2)
        KEYCOUNT[] += 1
        next_key1 = KEYCOUNT[] + 1
        next_key2 = KEYCOUNT[] + 2
    end

    # add two new entries for the bond
    s.pos[next_key1] = Vector{Int32}(undef, 5)
    s.pos[next_key1][1] = s.newPosi[1]
    s.pos[next_key1][2] = s.newPosi[2]
    s.pos[next_key1][3] = s.newPosi[3]
    s.pos[next_key1][4] = s.ki
    s.pos[next_key1][5] = next_key2

    s.pos[next_key2] = Vector{Int32}(undef, 5)
    s.pos[next_key2][1] = s.newPosj[1]
    s.pos[next_key2][2] = s.newPosj[2]
    s.pos[next_key2][3] = s.newPosj[3]
    s.pos[next_key2][4] = next_key1
    s.pos[next_key2][5] = s.kj

    # re-route the old pointers
    s.pos[s.ki][5] = next_key1
    s.pos[s.kj][4] = next_key2

    # Update keys vector
    push!(s.keys_vec, next_key1)
    push!(s.keys_vec, next_key2)

    KEYCOUNT[] += 2

    return s
end

export kinkBond!
function kinkBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi == s.bprev) ⊻ (s.newPosj == s.bnext) begin
        "Move is not a kink"
    end

    # set angled site to 0, flat site to 1 (allocation-free: no splatting or array construction)
    if s.newPosi == s.bprev 
        s.latS[s.bi[1], s.bi[2], s.bi[3]] = 0
        s.latS[s.newPosj[1], s.newPosj[2], s.newPosj[3]] = 1

        s.pos[s.ki] = Vector{Int32}(undef, 5)
        s.pos[s.ki][1] = s.newPosj[1]
        s.pos[s.ki][2] = s.newPosj[2]
        s.pos[s.ki][3] = s.newPosj[3]
        s.pos[s.ki][4] = s.kprev
        s.pos[s.ki][5] = s.kj

    elseif s.newPosj == s.bnext
        s.latS[s.newPosi[1], s.newPosi[2], s.newPosi[3]] = 1
        s.latS[s.bj[1], s.bj[2], s.bj[3]] = 0

        s.pos[s.kj] = Vector{Int32}(undef, 5)
        s.pos[s.kj][1] = s.newPosi[1]
        s.pos[s.kj][2] = s.newPosi[2]
        s.pos[s.kj][3] = s.newPosi[3]
        s.pos[s.kj][4] = s.ki
        s.pos[s.kj][5] = s.knext
    end

    return s
end

export removeBond!
function removeBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi == s.bprev) && (s.newPosj == s.bnext) begin
        "Move is not a bond removal"
    end
    
    # reset old lattice sites to 0
    s.latS[s.bi...] = 0
    s.latS[s.bj...] = 0

    # join ends together
    s.pos[s.kprev][5] = s.knext
    s.pos[s.knext][4] = s.kprev

    if s.pos[s.kprev][4] == 0
        s.head = s.kprev   # TEST
    end

    # Update keys vector - remove both keys
    filter!(k -> k != s.ki && k != s.kj, s.keys_vec)

    # remove keys associated with old bond
    delete!(s.pos, s.ki)
    delete!(s.pos, s.kj)

    return s
end

export identifyBond!
function identifyBond!(s::Species)

    @assert s.moveIdx >= 1 && s.moveIdx <= s.length begin
        "s.moveIdx must be between 1 and Species length (# of monomers)"
    end

    # keys_array = collect(keys(s.pos))  # TODO: this wastes allocs, fix it

    # pointers to the two beads (i and j) to be moved
    s.ki = s.keys_vec[s.moveIdx]
    s.kj = s.pos[s.ki][5]

    # previous and next beads along the chain
    s.kprev = s.pos[s.ki][4]
    s.knext = s.kj == 0 ? 0 : s.pos[s.kj][5]

    # bond edge coordinates
    s.bi = view(s.pos[s.ki], RANGE_1_3)
    s.bj = s.kj == 0 ? view(SENTINEL, RANGE_1_3) : view(s.pos[s.kj], RANGE_1_3)

    # coordinates of previous and next beads relative to bond
    s.bprev = s.kprev == 0 ? view(SENTINEL, RANGE_1_3) : view(s.pos[s.kprev], RANGE_1_3)
    s.bnext = s.knext == 0 ? view(SENTINEL, RANGE_1_3) : view(s.pos[s.knext], RANGE_1_3)

    for i in RANGE_1_3
        s.bondVec[i] = (s.bj[i] - s.bi[i]) % s.boxlength
    end

    return s
end

# CONTEXT: needs s.bondVec to be set
export moveTypes!
function moveTypes!(s::Species)

    # dimensions along which to move
    dims = findall(iszero, s.bondVec)

    @assert length(dims) == 2 begin
        "bondVec must be aligned with x, y, or z axis"
    end

    # build all possible moves each bond can make
    for (i, idx, sign) in zip(RANGE_1_4, COORD_MOVES, MOVE_ORIENTATIONS)
        
        for i in RANGE_1_3
            s.moveVec[i] = 0
        end

        # flip direction along each dimension in the orthogonal plane
        s.moveVec[dims[idx]] = sign

        # Compute bi + moveVec manually (no allocation)
        for j in RANGE_1_3
            s.tempVec[j] = s.bi[j] + s.moveVec[j]
        end
        boxmod!(s.tempVec, s.newPosi, Int32(1), s.boxlength)
        
        # Compute bj + moveVec manually (no allocation)
        for j in RANGE_1_3
            s.tempVec[j] = s.bj[j] + s.moveVec[j]
        end
        boxmod!(s.tempVec, s.newPosj, Int32(1), s.boxlength)
        
        
        # 1 -> bond addition, 2 -> kink, 3 -> bond removal
        s.moveTypes[i] = Int32(s.newPosi == s.bprev) + Int32(s.newPosj == s.bnext) + 1
    end

    for i in RANGE_1_3
        s.moveVec[i] = 0
    end

    return s
end

export makeMove!
function makeMove!(s::Species)

    @assert s.moveType == 1 || s.moveType == 2 || s.moveType == 3 begin
        "s.moveType must be either 1, 2, or 3"
    end

    # bond addition
    if s.moveType == 1
        addBond!(s)

    # kink
    elseif s.moveType == 2
        kinkBond!(s)

    # bond removal  
    elseif s.moveType == 3 
        removeBond!(s)
    end

    return s
end

# executes a move on polymer s inside lattice l at position s.moveIdx, of type
# moveType (add, remove or kink) sampled according to the metropolis algorithm
export movePoly!
function movePoly!(s::Species)
    
    @assert s.name[1] == 'p'
    moveExecuted = false

    # sample a bond at random along the chain
    s.moveIdx = rand(1:s.length)
    
    # find coordinates and keys of selected bond and adjacent beads
    identifyBond!(s)
    
    # don't move if last bead was selected
    if s.kj != 0
        
        # compute the possible move types around a bond
        moveTypes!(s)
        
        # sample a move according to proposition probabilities
        proposeMove!(s)
        
        # check for max length constraint and overlap
        if !overlap(s) && !(s.moveType == 1 && s.length == s.maxlength)
            
            # accept move with metropolis criterion
            if acceptMove(s)
                
                makeMove!(s)
                s.length = length(keys(s.pos))   # TODO: refactor for speed
                moveExecuted = true
            end
        end
    end

    return moveExecuted
end

# Executes a move on a species of bulk monomers diluted in the box
# by a 3D random walk on each particle
export moveBulkRw!
function moveBulkRw!(s::Species)
    
    @assert s.name[1] == 'b'
    moveExecuted = false

    # sample a bond at random along the chain
    s.moveIdx = rand(1:s.length)
    
    # turn off lattice site
    s.latS[s.pos[s.moveIdx][1], s.pos[s.moveIdx][2], s.pos[s.moveIdx][3]] = 0

    # sample an x/y/z direction and orientation at random
    s.moveType = rand(RANGE_1_3)
    s.pos[s.moveIdx][s.moveType] += 2 * Int(rand(Bernoulli(0.5))) - 1
    
    # clip to within bounding box
    if s.pos[s.moveIdx][s.moveType] > s.boxlength
        s.pos[s.moveIdx][s.moveType] -= s.boxlength
    elseif s.pos[s.moveIdx][s.moveType] < 0
        s.pos[s.moveIdx][s.moveType] += s.boxlength
    end

    # update new lattice site
    s.latS[s.pos[s.moveIdx][1], s.pos[s.moveIdx][2], s.pos[s.moveIdx][3]] = 1

    moveExecuted = true

    return moveExecuted
end
