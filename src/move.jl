
############# Polymer move functions #############

export overlap
function overlap(p::Poly)
    return ((p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] == 1 && 
             !coords_equal(p.bprev, p.newPosi)) ||
            (p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] == 1 && 
             !coords_equal(p.bnext, p.newPosj)))
end

export addBond!
function addBond!(p::Poly)

    # TODO: handle corner cases
    @assert (p.newPosi != p.bprev) && (p.newPosj != p.bnext) begin
        "Move is not a bond addition"
    end

    # set new lattice sites to 1
    p.latS[p.newPosi...] = 1
    p.latS[p.newPosj...] = 1

    # Find next available keys
    next_key1 = KEYCOUNT[] + 1
    next_key2 = KEYCOUNT[] + 2
    while haskey(p.pos, next_key1) || haskey(p.pos, next_key2)
        KEYCOUNT[] += 1
        next_key1 = KEYCOUNT[] + 1
        next_key2 = KEYCOUNT[] + 2
    end

    # add two new entries for the bond
    p.pos[next_key1] = Vector{Int32}(undef, 5)
    p.pos[next_key1][1] = p.newPosi[1]
    p.pos[next_key1][2] = p.newPosi[2]
    p.pos[next_key1][3] = p.newPosi[3]
    p.pos[next_key1][4] = p.ki
    p.pos[next_key1][5] = next_key2

    p.pos[next_key2] = Vector{Int32}(undef, 5)
    p.pos[next_key2][1] = p.newPosj[1]
    p.pos[next_key2][2] = p.newPosj[2]
    p.pos[next_key2][3] = p.newPosj[3]
    p.pos[next_key2][4] = next_key1
    p.pos[next_key2][5] = p.kj

    # re-route the old pointers
    p.pos[p.ki][5] = next_key1
    p.pos[p.kj][4] = next_key2

    # Update keys vector
    push!(p.keys_vec, next_key1)
    push!(p.keys_vec, next_key2)

    KEYCOUNT[] += 2

    return p
end

export kinkBond!
function kinkBond!(p::Poly)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) ⊻ (p.newPosj == p.bnext) begin
        "Move is not a kink"
    end

    # set angled site to 0, flat site to 1 (allocation-free: no splatting or array construction)
    if p.newPosi == p.bprev 
        p.latS[p.bi[1], p.bi[2], p.bi[3]] = 0
        p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] = 1

        p.pos[p.ki] = Vector{Int32}(undef, 5)
        p.pos[p.ki][1] = p.newPosj[1]
        p.pos[p.ki][2] = p.newPosj[2]
        p.pos[p.ki][3] = p.newPosj[3]
        p.pos[p.ki][4] = p.kprev
        p.pos[p.ki][5] = p.kj

    elseif p.newPosj == p.bnext
        p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] = 1
        p.latS[p.bj[1], p.bj[2], p.bj[3]] = 0

        p.pos[p.kj] = Vector{Int32}(undef, 5)
        p.pos[p.kj][1] = p.newPosi[1]
        p.pos[p.kj][2] = p.newPosi[2]
        p.pos[p.kj][3] = p.newPosi[3]
        p.pos[p.kj][4] = p.ki
        p.pos[p.kj][5] = p.knext
    end

    return p
end

export removeBond!
function removeBond!(p::Poly)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) && (p.newPosj == p.bnext) begin
        "Move is not a bond removal"
    end
    
    # reset old lattice sites to 0
    p.latS[p.bi...] = 0
    p.latS[p.bj...] = 0

    # join ends together
    p.pos[p.kprev][5] = p.knext
    p.pos[p.knext][4] = p.kprev

    if p.pos[p.kprev][4] == 0
        p.head = p.kprev   # TEST
    end

    # Update keys vector - remove both keys
    filter!(k -> k != p.ki && k != p.kj, p.keys_vec)

    # remove keys associated with old bond
    delete!(p.pos, p.ki)
    delete!(p.pos, p.kj)

    return p
end

# Finds the bead coordinates, neighbording beads, and vector of the bond 
# to which the bead selected by p.moveIdx belongs
export identifyBond!
function identifyBond!(p::Poly)

    @assert p.moveIdx >= 1 && p.moveIdx <= p.length begin
        "p.moveIdx must be between 1 and Species length (# of monomers)"
    end

    # keys_array = collect(keys(p.pos))  # TODO: this wastes allocs, fix it

    # pointers to the two beads (i and j) to be moved
    p.ki = p.keys_vec[p.moveIdx]
    p.kj = p.pos[p.ki][5]

    # previous and next beads along the chain
    p.kprev = p.pos[p.ki][4]
    p.knext = p.kj == 0 ? 0 : p.pos[p.kj][5]

    # bond edge coordinates
    p.bi = view(p.pos[p.ki], range(1, 3))
    p.bj = p.kj == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.kj], range(1, 3))

    # coordinates of previous and next beads relative to bond
    p.bprev = p.kprev == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.kprev], range(1, 3))
    p.bnext = p.knext == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.knext], range(1, 3))

    for i in range(1, 3)
        p.bondVec[i] = (p.bj[i] - p.bi[i]) % p.boxlength
    end

    return p
end

# Computes all possible move types for the selected bond
# Warning: needs p.bondVec to be set
export moveTypes!
function moveTypes!(p::Poly)

    # dimensions along which to move
    dims = findall(iszero, p.bondVec)

    @assert length(dims) == 2 begin
        "bondVec must be aligned with x, y, or z axis"
    end

    # build all possible moves each bond can make
    for (i, idx, sign) in zip(range(1, 4), COORD_MOVES, MOVE_ORIENTATIONS)
        
        for i in range(1, 3)
            p.moveVec[i] = 0
        end

        # flip direction along each dimension in the orthogonal plane
        p.moveVec[dims[idx]] = sign

        # Compute bi + moveVec
        for j in range(1, 3)
            p.tempVec[j] = p.bi[j] + p.moveVec[j]
        end
        boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)
        
        # Compute bj + moveVec
        for j in range(1, 3)
            p.tempVec[j] = p.bj[j] + p.moveVec[j]
        end
        boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)
        
        
        # 1 -> bond addition, 2 -> kink, 3 -> bond removal
        p.moveTypes[i] = Int32(p.newPosi == p.bprev) + Int32(p.newPosj == p.bnext) + 1
    end

    for i in range(1, 3)
        p.moveVec[i] = 0
    end

    return p
end

# Calls the appropriate operation (add, kink or remove) to execute a move
export makeMove!
function makeMove!(p::Poly)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    # bond addition
    if p.moveType == 1
        addBond!(p)

    # kink
    elseif p.moveType == 2
        kinkBond!(p)

    # bond removal  
    elseif p.moveType == 3 
        removeBond!(p)
    end

    return p
end

# executes a move on polymer s inside lattice l at position p.moveIdx, of type
# moveType (add, remove or kink) sampled according to the metropolis algorithm
export move!
function move!(p::Poly)

    # sample a bond at random along the chain
    p.moveIdx = rand(1:p.length)
    
    # find coordinates and keys of selected bond and adjacent beads
    identifyBond!(p)
    
    # don't move if last bead was selected
    if p.kj != 0
        
        # compute the possible move types around a bond
        moveTypes!(p)
        
        # sample a move according to proposition probabilities
        proposeMove!(p)
        
        # check for max length constraint and overlap
        if !overlap(p) && !(p.moveType == 1 && p.length == p.maxlength)
            
            # accept move with metropolis criterion
            if acceptMove(p)
                
                makeMove!(p)
                p.length = length(keys(p.pos))   # TODO: refactor for speed
                return true
            end
        end
    end

    return false
end


############# Bulk move functions #############

export overlap
function overlap(b::AbstractBulk)
    return b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] > 0
end

export identifyBond!
function identifyBond!(b::Bulk)
    b.bi = view(b.pos[b.moveIdx], range(1, 3))

    return b
end

# Translates a particle to a new position set by b.newPos
export makeMove!
function makeMove!(b::Bulk)
    
    # update lattice sites
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = 0
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # update current position
    for i in range(1, 3)
        b.pos[b.moveIdx][i] = b.newPos[i]
    end

    @assert b.latS[b.pos[b.moveIdx][1], 
                   b.pos[b.moveIdx][2], 
                   b.pos[b.moveIdx][3]] == b.moveIdx begin
        "Invalid index-position correspondence"
    end

    return b
end

# Executes a move on a species of bulk monomers diluted in the box
export move!
function move!(b::Bulk)

    # sample a particle at random
    b.moveIdx = rand(1:b.length)

    # find coordinates and keys of selected bond
    identifyBond!(b)
    
    # sample a move uniformly in each direction
    proposeMove!(b)

    if !overlap(b)
        if acceptMove(b)
            makeMove!(b)
            return true
        end
    end

    return false
end


############# FastBulk move functions #############

export identifyBond!
function identifyBond!(b::FastBulk)
    b.bi = view(b.pos[b.moveIdx], range(1, 3))

    return b
end

# Translates a particle to a new position set by b.newPos
export makeMove!
function makeMove!(b::FastBulk)
    
    # delete old lattice site
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = 0

    # find current nns, make them movable
    nnIdxs!(b, b.bi)
    
    for nnIdx in b.nnIdxs
        if nnIdx > 0
            push!(b.movableIdxs, nnIdx)
        end
    end

    # set new lattice site
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # find if new position's nns are no longer movable
    nnIdxs!(b, b.newPos)

    for nnIdx in b.nnIdxs
        if nnIdx > 0
            if nn(b, b.pos[nnIdx][range(1, 3)]) == 6
                pop!(b.movableIdxs, nnIdx)
            end
        end
    end

    # add new site to movableIdxs (most be movable to preserve ergodicity)
    push!(b.movableIdxs, b.moveIdx)
    
    # update current position and new lattice site
    for i in range(1, 3)
        b.pos[b.moveIdx][i] = b.newPos[i]
    end

    @assert b.latS[b.pos[b.moveIdx][1], 
                   b.pos[b.moveIdx][2], 
                   b.pos[b.moveIdx][3]] == b.moveIdx begin
        "Invalid index-position correspondence"
    end

    return b
end

# Executes a move on a species of bulk monomers diluted in the box
export move!
function move!(b::FastBulk)

    # sample a particle at random among movable particles
    b.moveIdx = rand(b.movableIdxs)

    # find coordinates and keys of selected bond
    identifyBond!(b)
    
    # sample a move uniformly in each direction
    proposeMove!(b)

    @assert !overlap(b)

    if acceptMove(b)
        makeMove!(b)
        return true
    end

    return false
end
