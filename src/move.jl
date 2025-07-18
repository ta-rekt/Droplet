# executes a move on p at position movePos, in a direction
# specified by moveDir, which specifies the axes along which the move happens
export makeMove!
function makeMove!(p::Poly, movePos::Int, moveDir::Vector{T}) where T <: AbstractFloat
    @assert length(moveDir) == 3 "moveDir must be a vector of length 3"
    @assert movePos >= 1 || movePos <= p.length "movePos must be between 1 and polymer length"

    # bond to move is between beads i and j
    bi = view(p.pos, movePos, :)
    bj = view(p.pos, (movePos + 1) % p.boxlength, :)
    bond_vec = [(bj[1] - bi[1]) % p.boxlength, 
                (bj[2] - bi[2]) % p.boxlength, 
                (bj[3] - bi[3]) % p.boxlength]

    # previous and next bead relative to selected bond
    ni = [-1, -1, -1]
    nj = [-1, -1, -1]

    if movePos != 1
        ni = view(p.pos, movePos - 1, :)
    end

    if movePos != p.length - 1
        nj = view(p.pos, movePos + 2, :)
    end

    # 

end

function proposeMove!(p::Poly)

end

function acceptMove!(p::Poly)

end