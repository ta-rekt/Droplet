# samples moves according to probabilities set by detailed balance
# writes the sampled direction vector in s.moveVec
export proposeMove!
function proposeMove!(s::Species)

    # compute how many of each move type there are
    nTypes = zeros(Int32, 3)

    @assert !isnothing(s.moveTypes) begin
        "moveTypes! has not been called yet"
    end

    for i in 1:4
        nTypes[s.moveTypes[i]] += 1
    end

    # set Pkink to zero if no kinks are possible
    Pkink = (nTypes[2] == 0) ? 0 : s.moveProbs[2]

    N = s.length
    dims = findall(iszero, s.bondVec)
    
    s.moveProbs[1] = (1 - Pkink) * (N + 2) / (2*N + 2)
    s.moveProbs[3] = (1 - Pkink) * N / (2*N + 2)

    # sample move types according to proposed probabilities
    r = rand(Uniform(0,1))
    p = 0.0

    for (i, idx, sign) in zip(1:4, [1,1,2,2], [+1,-1,+1,-1])
        tp = s.moveTypes[i]
        p += s.moveProbs[tp] / nTypes[tp]

        if r < p
            s.moveVec[dims[idx]] = sign
            s.moveType = tp
            
            return s
        end
    end

    @assert dot(s.bondVec, s.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end

    @assert p == 1 begin
        "move probs don't sum to one"
    end

    return s
end

# accepts of rejects according to the metropolis criterion
export acceptMove!
function acceptMove!(s::Species)
    
    return true
end