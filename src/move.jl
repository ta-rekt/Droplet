# executes a move on species s inside lattice l at position movePos
# in a direction moveVec, the vector along which the move happens
export makeMove!
function makeMove!(l::Lattice, 
                   s::Species, 
                   movePos::Int, 
                   moveVec::Vector{T}) where T <: Int

    @assert length(moveVec) == 3 
        "moveVec must be a vector of length 3"
    @assert movePos >= 1 || movePos <= s.length 
        "movePos must be between 1 and Species length (# of monomers)"

    # lattice sites of species s
    latS = view(l.lat, l.speciesIndex[s.id], :, :, :)

    # bond to move is between beads i and j
    bi = view(s.pos, movePos, :)
    bj = view(s.pos, (movePos + 1) % s.boxlength, :)
    s.bondVec = [(bj[1] - bi[1]) % s.boxlength, 
                 (bj[2] - bi[2]) % s.boxlength, 
                 (bj[3] - bi[3]) % s.boxlength]

    # check move is orthogonal to bond direction
    @assert dot(s.bondVec, moveVec) == 0 
        "Move is not orthogonal to bond direction"

    # compute translated position of selected bond
    for k in 1:3
        s.auxPos1[k] = bi[k] + moveVec[k]
        s.auxPos2[k] = bj[k] + moveVec[k]
    end
    
    # check for overlap with another bead from the same species
    # WARNING: different species can still overlap
    if latS[CartesianIndex((s.auxPos1[1], s.auxPos1[2], s.auxPos1[3]))] != 1
        if latS[CartesianIndex((s.auxPos2[1], s.auxPos2[2], s.auxPos2[3]))] != 1

            # reset old lattice sites to 0
            # latS[CartesianIndex((bi[1], bi[2], bi[3]))] = 0
            # latS[CartesianIndex((bj[1], bj[2], bj[3]))] = 0
            
            # flip new lattice sites to 1
            latS[CartesianIndex((s.auxPos1[1], s.auxPos1[2], s.auxPos1[3]))] = 1
            latS[CartesianIndex((s.auxPos2[1], s.auxPos2[2], s.auxPos2[3]))] = 1
            
            # update bead position 
            for k in 1:3
                bi[k] = s.auxPos1[k]
                bj[k] = s.auxPos2[k]
            end
        end
    end

    return l, s
end

function proposeMove!(s::Species)

end

function acceptMove!(s::Species)

end

# # previous and next bead relatSve to selected bond
# if movePos != 1
#     ni = view(s.pos, movePos - 1, :)
# else
#     ni = [-1, -1, -1]
# end

# if movePos != s.length - 1
#     nj = view(s.pos, movePos + 2, :)
# else
#     nj = [-1, -1, -1]
# end