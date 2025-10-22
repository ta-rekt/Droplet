# samples moves according to probabilities set by detailed balance
# writes the sampled direction vector in s.moveVec
# CONTEXT: needs s.moveTypes to be set
export proposeMove!
function proposeMove!(s::Species)
    
    @assert !isnothing(s.moveTypes) begin
        "moveTypes! has not been called yet"
    end

    # compute how many of each move type there are
    s.nTypes = zeros(Int32, 3)

    for i in 1:4
        s.nTypes[s.moveTypes[i]] += 1
    end
    
    N = s.length
    
    # set Pkink to zero if no kinks are possible
    Pkink = s.nTypes[2] == 0 ? 0 : s.moveTypeProbs[2]
    
    # Padd
    s.moveTypeProbs[1] = s.nTypes[3] == 0 ? (1 - Pkink) : (1 - Pkink) * (N + 2) / (2*N + 2)
    
    # Premove
    s.moveTypeProbs[3] = s.nTypes[3] == 0 ? 0 : (1 - Pkink) * N / (2*N + 2)

    @assert s.moveTypeProbs[1] + s.moveTypeProbs[3] + Pkink == 1 begin
        "Proposition probabilities don't sum to 1"
    end

    dims = findall(iszero, s.bondVec)

    # sample move types according to proposed probabilities
    r = rand(Uniform(0,1))
    p = 0.0

    for (i, idx, sign) in zip(1:4, [1,1,2,2], [+1,-1,+1,-1])
        
        s.moveVec = [0,0,0]
        typ = s.moveTypes[i]
        p += s.moveTypeProbs[typ] / s.nTypes[typ]

        if r < p
            s.moveVec[dims[idx]] = sign
            s.moveType = typ
            break
        end
    end

    @assert dot(s.bondVec, s.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end

    # translated position of selected bond
    # s.newPosi = s.bi .+ s.moveVec
    # s.newPosj = s.bj .+ s.moveVec

    boxmod!(s.bi .+ s.moveVec, s.newPosi, Int32(1), s.boxlength)
    boxmod!(s.bj .+ s.moveVec, s.newPosj, Int32(1), s.boxlength)

    return s
end

# accepts or rejects according to the metropolis criterion
export acceptMove
function acceptMove(s::Species)

    deltaHamPoly!(s)
    
    p = min(1, exp(-s.deltaH))
    r = rand(Uniform(0,1))

    if r < p
        return true
    else
        return false
    end
end