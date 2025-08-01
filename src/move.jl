export removeBond!
function removeBond!(latS::S, s::Species, ki::Int32, kj::Int32, kprev::Int32, knext::Int32,
                     bi::T, bj::T, bprev::T, bnext::T) where {T <: AbstractArray{Int32}, S <: AbstractArray{Int32}}

    # TODO: handle corner cases
    @assert (s.newPosi == bprev) && (s.newPosj == bnext) begin
        "Move is not a bond removal"
    end
    
    # reset old lattice sites to 0
    latS[bi...] = 0
    latS[bj...] = 0

    # join ends together
    s.pos[kprev][5] = knext
    s.pos[knext][4] = kprev

    # remove keys associated with old bond
    delete!(s.pos, ki)
    delete!(s.pos, kj)

    s.length = s.length - 2
end

export addBond!
function addBond!(latS::S, s::Species, ki::Int32, kj::Int32, kprev::Int32, knext::Int32,
                  bi::T, bj::T, bprev::T, bnext::T) where {T <: AbstractArray{Int32}, S <: AbstractArray{Int32}}

    # TODO: handle corner cases
    @assert (s.newPosi != bprev) && (s.newPosj != bnext) begin
        "Move is not a bond addition"
    end

    # set new lattice sites to 1
    latS[s.newPosi...] = 1
    latS[s.newPosj...] = 1

    # add two new entries for the bond
    s.pos[s.length+1] = [s.newPosi..., ki, s.length+2]
    s.pos[s.length+2] = [s.newPosj..., s.length+1, kj]

    # re-route the old pointers
    s.pos[ki][5] = s.length + 1
    s.pos[kj][4] = s.length + 2

    s.length = s.length + 2
end

export kinkBond!
function kinkBond!(latS::S, s::Species, ki::Int32, kj::Int32, kprev::Int32, knext::Int32,
                   bi::T, bj::T, bprev::T, bnext::T) where {T <: AbstractArray{Int32}, S <: AbstractArray{Int32}}

    # TODO: handle corner cases
    @assert (s.newPosi == bprev) ⊻ (s.newPosj == bnext) begin
        "Move is not a kink"
    end

    # set angled site to 0, flat site to 1
    if s.newPosi == bprev 
        latS[bi...] = 0
        latS[s.newPosj...] = 1

        s.pos[ki] = [s.newPosj..., kprev, kj]

    elseif s.newPosj == bnext
        latS[s.newPosi...] = 1
        latS[bj...] = 0

        s.pos[kj] = [s.newPosi..., ki, knext]
    end
end

export moveTypes!
function moveTypes!(latS::S, s::Species, ki::Int32, kj::Int32, kprev::Int32, knext::Int32,
                  bi::T, bj::T, bprev::T, bnext::T) where {T <: AbstractArray{Int32}, S <: AbstractArray{Int32}}

    

end

# executes a move on species s inside lattice l at position s.moveIdx,
# of type moveType (add, remove or kink) in a direction s.moveVec, 
# the vector along which the move happens
export makeMove!
function makeMove!(l::Lattice, s::Species)

    @assert s.moveIdx >= 1 && s.moveIdx <= s.length - 1 begin
        "s.moveIdx must be between 1 and Species length - 1 (# of monomers)"
    end

    moveExecuted = false

    # lattice sites of species s
    latS = view(l.lat, l.speciesDir[s.id], :, :, :)

    # pointers to the two beads (i and j) to be moved, 
    ki = s.pos.keys[s.moveIdx]
    kj = s.pos[ki][5]

    # previous and next beads along the chain
    kprev = s.pos[ki][4]
    knext = s.pos[kj][5]

    # bond edge coordinates
    bi = view(s.pos[ki], 1:3)
    bj = view(s.pos[kj], 1:3)

    # coordinates of previous and next beads relative to bond
    bprev = kprev == 0 ? view(Int32[-1,-1,-1], 1:3) : view(s.pos[kprev], 1:3)
    bnext = knext == 0 ? view(Int32[-1,-1,-1], 1:3) : view(s.pos[knext], 1:3)

    s.bondVec = [(bj[1] - bi[1]) % s.boxlength, 
                 (bj[2] - bi[2]) % s.boxlength, 
                 (bj[3] - bi[3]) % s.boxlength]
    
    ###########

    # sample a move, its type and direction are in s's fields
    proposeMove!(s)

    @assert !isnothing(s.moveVec)
    @assert length(s.moveVec) == 3 begin
        "s.moveVec must be a vector of length 3"
    end
    @assert dot(s.bondVec, s.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end
    @assert s.moveType == 1 || s.moveType == 2 || s.moveType == 3 begin
        "s.moveType must be either 1, 2, or 3"
    end

    # translated position of selected bond
    s.newPosi = bi .+ s.moveVec
    s.newPosj = bj .+ s.moveVec

    # check for max length constraint and overlap
    if !((s.moveType == 2 && s.length == s.maxlength) ||
        (latS[s.newPosi...] == 1 && bprev != s.newPosi) ||
        (latS[s.newPosj...] == 1 && bnext != s.newPosj))

        # bond removal
        if s.moveType == 1  
            removeBond!(latS, s, ki, kj, kprev, knext, bi, bj, bprev, bnext)

        # bond addition
        elseif s.moveType == 2  
            addBond!(latS, s, ki, kj, kprev, knext, bi, bj, bprev, bnext)

        # kink
        elseif s.moveType == 3 
            kinkBond!(latS, s, ki, kj, kprev, knext, bi, bj, bprev, bnext)
        end

        moveExecuted = acceptMove!(s)
    end

    return l, s, moveExecuted
end