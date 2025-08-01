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

# executes a move on species s inside lattice l at position s.moveIdx,
# of type moveType (add, remove or kink) in a direction s.moveVec, 
# the vector along which the move happens
export makeMove!
function makeMove!(l::Lattice, s::Species)

    @assert !isnothing(s.moveVec)
    @assert length(s.moveVec) == 3 begin
        "s.moveVec must be a vector of length 3"
    end
    @assert s.moveIdx >= 1 && s.moveIdx <= s.length - 1 begin
        "s.moveIdx must be between 1 and Species length - 1 (# of monomers)"
    end
    @assert s.moveType == 1 || s.moveType == 2 || s.moveType == 3 begin
        "s.moveType must be either 1, 2, or 3"
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

    @assert dot(s.bondVec, s.moveVec) == 0 begin
        "Move is not orthogonal to bond direction"
    end

    # translated position of selected bond
    s.newPosi = bi .+ s.moveVec
    s.newPosj = bj .+ s.moveVec

    # check for length constraint and overlap
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

        moveExecuted = true
    end

    return l, s, moveExecuted
end


# selects a bond randomly, proposes a move following detailed balance
# checks for overlap and length constraints, and accepts with metropolis prob.
function sampleMove!(l::Lattice, s::Species)
    # select a bond randomly
    s.moveIdx = rand(1:(s.length-1))

    # propose a move following detailed balance

    # check for overlap and length constraints

    # accept with metropolis probability

    # return l, s
end

# moves a random bead according to the metropolis rule
export move!
function move!(l::Lattice, s::Species)

    sampleMove!(l::Lattice, s::Species)
    makeMove!(l::Lattice, s::Species)

    return l, s
end

# old code for array implementation of s.pos

    # bi = view(s.pos, s.moveIdx, :)
    # bj = view(s.pos, (s.moveIdx + 1) % s.boxlength, :)
    # bprev = s.moveIdx == 1 ? [-1, -1, -1] : view(s.pos, s.moveIdx - 1, :)
    # bnext = s.moveIdx == s.length - 1 ? [-1, -1, -1] : view(s.pos, s.moveIdx + 2, :)

    ### 

    # # reset old lattice sites to 0
    # latS[bi...] = 0
    # latS[bj...] = 0

    # # move elements of pos array back by 2 indices
    # s.pos[s.moveIdx:end-2, :] = s.pos[s.moveIdx+2:end, :]
    # s.pos[end-1, :] = zeros(Int, 3)
    # s.pos[end, :] = zeros(Int, 3)

    ## 

    # s.length = s.moveIdx == 1 || s.moveIdx == (s.length - 1) ? s.length -1 : s.length - 2
    # # set new lattice sites to 1
    # latS[s.newPosi...] = 1
    # latS[s.newPosj...] = 1

    # # move elements of pos forward by 2 indices
    # for i in (s.maxlength-2):-1:s.moveIdx
    #     rowi = view(s.pos, i, :)
    #     rowj = view(s.pos, i+2, :)
    #     rowj = rowi
    # end
    
    ## 

    # # update bead position
    # bi = s.newPosi
    # bj = s.newPosj

    # s.length = s.moveIdx == 1 || s.moveIdx == (s.length - 1) ? s.length + 1 : s.length + 2
    # # set angled site to 0, flat site to 1
    # if s.newPosi == bprev
    #     latS[bi...] = 0
    #     latS[s.newPosj...] = 1

    # elseif s.newPosj == bnext
    #     latS[s.newPosi...] = 1
    #     latS[bj...] = 0
    # end