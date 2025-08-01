# samples moves according to probabilities set by detailed balance
export proposeMove!
function proposeMove!(s::Species)

    Pkink = s.moveProbs[3]
    N = s.length

    @assert Pkink > 0 begin
        "Pkink must be nonzero"
    end
    
    s.moveProbs[1] = (1 - Pkink) * N / (2*N + 2)
    s.moveProbs[2] = (1 - Pkink) * (N + 2) / (2*N + 2)

    # sample a move

end

# accepts of rejects according to the metropolis criterion
export acceptMove!
function acceptMove!(s::Species)

end