############# Generic metropolis functions #############
# MARK: -Generic

# accepts or rejects according to the metropolis criterion
export acceptMove
function acceptMove(s::AbstractSpecies)

    deltaHam!(s)
    prob = min(1, s.hastingsFactor * exp(-s.deltaH))

    return rand(Uniform(0,1)) < prob
end


############# Polymer metropolis functions #############
# MARK: -Poly

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
        p.moveTypes[i] = (Int32(p.newPosi == p.bprev) 
                        + Int32(p.newPosj == p.bnext) + 1)
    end

    for i in range(1, 3)
        p.moveVec[i] = 0
    end

    return p
end

# samples moves according to probabilities set by detailed balance
# writes the sampled direction vector in p.moveVec
export proposeMove!
function proposeMove!(p::Poly)
    
    # compute the possible move types around a bond
    moveTypes!(p)

    # compute how many of each move type there are
    for j in range(1, 3)
        p.nTypes[j] = 0
    end

    for i in range(1, 4)
        p.nTypes[p.moveTypes[i]] += 1
    end
    
    N = p.length
    
    # set Pkink to zero if no kinks are possible
    Pkink = p.nTypes[2] == 0 ? 0 : p.moveTypeProbs[2]
    
    # Padd
    p.moveTypeProbs[1] = p.nTypes[3] == 0 ? (1 - Pkink) : (1 - Pkink) * (N + 2) / (2*N + 2)
    
    # Premove
    p.moveTypeProbs[3] = p.nTypes[3] == 0 ? 0 : (1 - Pkink) * N / (2*N + 2)

    @assert p.moveTypeProbs[1] + p.moveTypeProbs[3] + Pkink == 1 begin
        "Proposition probabilities don't sum to 1"
    end

    dims = findall(iszero, p.bondVec)

    # sample move types according to proposed probabilities
    r = rand(Uniform(0,1))
    prob = 0.0

    for (i, idx, sign) in zip(range(1, 4), COORD_MOVES, MOVE_ORIENTATIONS)
        
        for i in range(1, 3)
            p.moveVec[i] = 0
        end
        
        typ = p.moveTypes[i]
        prob += p.moveTypeProbs[typ] / p.nTypes[typ]

        if r < prob
            p.moveVec[dims[idx]] = sign
            p.moveType = typ
            break
        end
    end

    @assert dot(p.bondVec, p.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end

    # new position is old position translated in proposed direction
    for j in range(1, 3)
        p.tempVec[j] = p.bi[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)
    
    for j in range(1, 3)
        p.tempVec[j] = p.bj[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)

    return p
end


############# Bulk metropolis functions #############
# MARK: -Bulk

export proposeMove!
function proposeMove!(b::Bulk)

    for i in range(1, 3)
        b.tempVec[i] = b.pos[b.moveIdx][i]
    end

    # sample an x/y/z direction and orientation at random
    b.moveType = rand(range(1, 3))
    b.tempVec[b.moveType] += 2 * Int(rand(Bernoulli(0.5))) - 1
    
    # clip to within bounding box
    boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)

    return b
end


############# FastBulk metropolis functions #############
# MARK: -FastBulk

export proposeMove!
function proposeMove!(b::FastBulk)

    nnIdxs!(b, b.bi)

    @assert !isempty(b.freeNNs) begin
        "No free nearest neighbor sites: proposing an impossible move"
    end

    for i in range(1, 3)
        b.tempVec[i] = b.pos[b.moveIdx][i]
    end

    # propose a move uniformly among the free sites
    b.moveType = rand(b.freeNNs)
    axis = div(b.moveType + 1, 2)
    dir = - mod(b.moveType + 1, 2) * 2 + 1
    b.tempVec[axis] += dir

    # clip to within bounding box
    boxmod!(b.tempVec, b.newPos, Int32(1), b.boxlength)
    
    hastingsFactor!(b)

    return b
end

# hastings factor. two terms: 
# 1. |S(x)| / |S(x')|, where S(x) is the movable set of current config. x
# 2. |M(x)| / |M(x')|, where M(x) is the number of free nns around b.bi
export hastingsFactor!
function hastingsFactor!(b::FastBulk)

    nmovable = length(b.movableIdxs)
    @assert nmovable > 0 "No movable beads"

    # (# of nnIdxs(bi) with nn = 6) - (# of nnIdxs(newPos) with nn = 5)
    nmovable_new = nmovable

    for key in b.nnIdxs
        if key != 0 && nn(b, b.pos[key][range(1,3)]) == 6
            nmovable_new += 1
        end
    end

    nnIdxs!(b, b.newPos)

    for key in b.nnIdxs
        if key != 0 && key != b.moveIdx && nn(b, b.pos[key][range(1,3)]) == 5
            nmovable_new -= 1
        end
    end

    nmoves = 6 - nn(b, b.bi)
    nmoves_new = 6 - (nn(b, b.newPos) - 1)  # count out current position
    @assert nmoves > 0 "Selected bead is not movable"

    b.hastingsFactor = nmovable * nmoves / (nmovable_new * nmoves_new)

    return b
end