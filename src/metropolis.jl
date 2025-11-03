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

# accepts or rejects according to the metropolis criterion
export acceptMove
function acceptMove(p::Poly)

    deltaHam!(p)
    
    p = min(1, exp(-p.deltaH))
    r = rand(Uniform(0,1))

    if r < p
        return true
    else
        return false
    end
end


############# Bulk metropolis functions #############

