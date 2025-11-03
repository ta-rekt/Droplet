############# Utilities #############

# Number of nearest neighbors at a given position
export nn
function nn(s::Species, pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += s.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]]

        pos[i] -= 2
        boxmod!(pos, s.tempVec, Int32(1), s.boxlength)
        n += s.latS[s.tempVec[1], s.tempVec[2], s.tempVec[3]]

        pos[i] += 1
    end

    return n
end


############# Polymer Energy functions #############

# Hamiltonian of a polymer
export hamPoly!
function hamPoly!(s::Species)

    s.ham = 0
    nnn = 0                # number of nearest neighbors
    key = s.head           # initial bead
    pos = s.pos[key][range(1, 3)]  # initial posiion

    for i in 1:s.length
        nnn += nn(s, pos)
        key = s.pos[key][5]    # next bead

        if key != 0
            pos = s.pos[key][range(1, 3)]  # next position
        end
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    s.ham = - (s.chemPot + s.intPot) * s.length + s.intPot * nnn

    return s
end

# Energy increment of a proposed move. Requires proposeMove! to be called
# which sets moveType and newPosi/j. 
export deltaHamPoly!
function deltaHamPoly!(s::Species)

    @assert s.moveType == 1 || s.moveType == 2 || s.moveType == 3 begin
        "s.moveType must be either 1, 2, or 3"
    end

    s.deltaH = 0

    # bond addition
    if s.moveType == 1
        s.deltaH = (- 2 * (s.chemPot + s.intPot) 
                    + s.intPot * (2*(nn(s, s.newPosi) + nn(s, s.newPosj)) + 2))

    # kink
    elseif s.moveType == 2

        if s.newPosi == s.bprev 
            s.deltaH = 2 * s.intPot * (nn(s, s.newPosj) - nn(s, s.bi))
        elseif s.newPosj == s.bnext
            s.deltaH = 2 * s.intPot * (nn(s, s.newPosi) - nn(s, s.bj))
        else
            error("Move is not a kink")
        end

    # bond removal  
    elseif s.moveType == 3 
        s.deltaH = (2 * (s.chemPot + s.intPot)  
                    - s.intPot * (2*(nn(s, s.bi) + nn(s, s.bj)) - 2))
    end

    s.ham += s.deltaH

    return s
end


############# Bulk energy functions #############

# Hamiltonian of a dilute species in the bulk
export hamBulk1s!
function hamBulk1s!(s::Species)
    
    s.ham = 0
    nnn = 0    # number of nearest neighbors

    for i in 1:s.length
        pos = s.pos[i][range(1, 3)]
        nnn += nn(s, pos)
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    s.ham = - s.chemPot * s.length + s.intPot * nnn

    return s
end

# Energy increment for the bulk system
export deltaHamBulk1s!
function deltaHamBulk1s!(s::Species)

    s.deltaH = 2 * s.intPot * (nn(s, s.newPosi) - nn(s, s.bi))
    s.ham += s.deltaH

    return s
end