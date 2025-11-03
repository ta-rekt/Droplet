############# Generic metropolis functions #############

# accepts or rejects according to the metropolis criterion
export acceptMove
function acceptMove(s::AbstractSpecies)

    deltaHam!(s)
    
    prob = min(1, exp(-s.deltaH))
    r = rand(Uniform(0,1))

    if r < prob
        return true
    else
        return false
    end
end


############# Polymer metropolis functions #############

# samples moves according to probabilities set by detailed balance
# writes the sampled direction vector in p.moveVec
# CONTEXT: needs p.moveTypes to be set
export proposeMove!
function proposeMove!(p::Poly)
    
    @assert !isnothing(p.moveTypes) begin
        "moveTypes! has not been called yet"
    end

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

    # translated position of selected bond (no allocation)
    for j in range(1, 3)
        p.tempVec[j] = p.bi[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)
    
    # translated position of selected bond (no allocation)
    for j in range(1, 3)
        p.tempVec[j] = p.bj[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)

    return p
end


############# Bulk metropolis functions #############

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