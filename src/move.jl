export addBond!
function addBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi != s.bprev) && (s.newPosj != s.bnext) begin
        "Move is not a bond addition"
    end

    # set new lattice sites to 1
    s.latS[s.newPosi...] = 1
    s.latS[s.newPosj...] = 1

    # add two new entries for the bond
    s.pos[s.length+1] = [s.newPosi..., s.ki, s.length+2]
    s.pos[s.length+2] = [s.newPosj..., s.length+1, s.kj]

    # re-route the old pointers
    s.pos[s.ki][5] = s.length + 1
    s.pos[s.kj][4] = s.length + 2

    s.length = s.length + 2

    return s
end

export kinkBond!
function kinkBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi == s.bprev) ⊻ (s.newPosj == s.bnext) begin
        "Move is not a kink"
    end

    # set angled site to 0, flat site to 1
    if s.newPosi == s.bprev 
        s.latS[s.bi...] = 0
        s.latS[s.newPosj...] = 1

        s.pos[s.ki] = [s.newPosj..., s.kprev, s.kj]

    elseif s.newPosj == s.bnext
        s.latS[s.newPosi...] = 1
        s.latS[s.bj...] = 0

        s.pos[kj] = [s.newPosi..., s.ki, s.knext]
    end

    return s
end

export removeBond!
function removeBond!(s::Species)

    # TODO: handle corner cases
    @assert (s.newPosi == s.bprev) && (s.newPosj == s.bnext) begin
        "Move is not a bond removal"
    end
    
    # reset old lattice sites to 0
    s.latS[s.bi...] = 0
    s.latS[s.bj...] = 0

    # join ends together
    s.pos[s.kprev][5] = s.knext
    s.pos[s.knext][4] = s.kprev

    # remove keys associated with old bond
    delete!(s.pos, s.ki)
    delete!(s.pos, s.kj)

    s.length = s.length - 2

    return s
end

export moveTypes!
function moveTypes!(s::Species)

    # dimensions along which to move
    dims = findall(iszero, s.bondVec)
    @assert length(dims) == 2 begin
        "bondVec must be aligned with x, y, or z axis"
    end

    # build all possible moves each bond can make
    for (i, idx, sign) in zip(1:4, [1,1,2,2], [+1,-1,+1,-1])
        
        # flip direction along each dimension in the orthogonal plane
        s.moveVec[dims[idx]] = sign
        
        s.newPosi = s.bi .+ s.moveVec
        s.newPosj = s.bj .+ s.moveVec
        
        # 1 -> bond addition, 2 -> kink, 3 -> bond removal
        s.moveTypes[i] = Int32(s.newPosi == s.bprev) + Int32(s.newPosj == s.bnext) + 1
    end

    return s
end

export identifyBond!
function identifyBond!(s::Species)

    @assert s.moveIdx >= 1 && s.moveIdx <= s.length - 1 begin
        "s.moveIdx must be between 1 and Species length - 1 (# of monomers)"
    end

    # pointers to the two beads (i and j) to be moved, 
    s.ki = s.pos.keys[s.moveIdx]
    s.kj = s.pos[s.ki][5]

    # previous and next beads along the chain
    s.kprev = s.pos[s.ki][4]
    s.knext = s.pos[s.kj][5]

    # bond edge coordinates
    s.bi = view(s.pos[s.ki], 1:3)
    s.bj = view(s.pos[s.kj], 1:3)

    # coordinates of previous and next beads relative to bond
    s.bprev = s.kprev == 0 ? view(Int32[-1,-1,-1], 1:3) : view(s.pos[s.kprev], 1:3)
    s.bnext = s.knext == 0 ? view(Int32[-1,-1,-1], 1:3) : view(s.pos[s.knext], 1:3)

    s.bondVec = [(s.bj[1] - s.bi[1]) % s.boxlength, 
                 (s.bj[2] - s.bi[2]) % s.boxlength, 
                 (s.bj[3] - s.bi[3]) % s.boxlength]

    return s
end

export makeMove!
function makeMove!(s::Species)

    @assert s.moveType == 1 || s.moveType == 2 || s.moveType == 3 begin
        "s.moveType must be either 1, 2, or 3"
    end

    # bond addition
    if s.moveType == 1
        addBond!(s)

    # kink
    elseif s.moveType == 2
        kinkBond!(s)

    # bond removal  
    elseif s.moveType == 3 
        removeBond!(s)
    end

    return s
end

# executes a move on species s inside lattice l at position s.moveIdx,
# of type moveType (add, remove or kink) in a direction s.moveVec, 
# the vector along which the move happens
export move!
function move!(l::Lattice, s::Species)
    
    moveExecuted = false

    # give species s a view on its lattice sites
    s.latS = view(l.lat, l.speciesDir[s.id], :, :, :)

    # find coordinates and keys of selected bond and adjacent beads
    identifyBond!(s)
    
    # compute the possible move types around a bond
    moveTypes!(s)

    # sample a move according to detailed balance
    proposeMove!(s)

    # translated position of selected bond
    s.newPosi = s.bi .+ s.moveVec
    s.newPosj = s.bj .+ s.moveVec

    # check for max length constraint and overlap
    if !((s.moveType == 2 && s.length == s.maxlength) ||
        (s.latS[s.newPosi...] == 1 && s.bprev != s.newPosi) ||
        (s.latS[s.newPosj...] == 1 && s.bnext != s.newPosj))

        # accept move with metropolis criterion
        if acceptMove!(s)
            
            makeMove!(s)
            moveExecuted = true
        end
    end

    return moveExecuted
end