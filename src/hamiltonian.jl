# Hamiltonian of a polymer
export hamPoly!
function hamPoly!(s::Species)

    s.ham = 0
    nnn = 0                # number of nearest neighbors
    key = s.head           # initial bead
    pos = s.pos[key][1:3]  # initial posiion

    for i in 1:s.length
        nnn += nn(s, pos)
        key = s.pos[key][5]    # next bead

        if key != 0
            pos = s.pos[key][1:3]  # next position
        end
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    # @show nnn

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
    # @show nn(s, s.newPosi)
    # @show nn(s, s.newPosj)

    # bond addition
    if s.moveType == 1
        # add 4 = 2 (canceled by 1st term) + 2 (move hasn't been made yet)
        s.deltaH = (- 2 * (s.chemPot + s.intPot) 
                    + s.intPot * (4 + nn(s, s.newPosi) + nn(s, s.newPosj)))

    # kink
    elseif s.moveType == 2

        if s.newPosi == s.bprev 
            s.deltaH = s.intPot * (nn(s, s.newPosj) - nn(s, s.bi))
        elseif s.newPosj == s.bnext
            s.deltaH = s.intPot * (nn(s, s.newPosi) - nn(s, s.bj))
        else
            error("Move is not a kink")
        end

    # bond removal  
    elseif s.moveType == 3 
        s.deltaH = (2 * (s.chemPot + s.intPot)  
                    - s.intPot * (2 + nn(s, s.bi) + nn(s, s.bj)))
    end

    s.ham += s.deltaH

    return s
end

# Number of nearest neighbors at a given position
export nn
function nn(s::Species, pos::Union{SubArray{Int32}, Vector{Int32}})
    n = 0
    adjPos = zeros(Int32, 3)

    for i in 1:3
        pos[i] += 1
        boxmod!(pos, adjPos, Int32(1), s.boxlength)
        n += s.latS[adjPos...]
        # @show adjPos
        # @show i, pos, n

        pos[i] -= 2
        boxmod!(pos, adjPos, Int32(1), s.boxlength)
        n += s.latS[adjPos...]
        # @show adjPos
        # @show i, pos, n

        pos[i] += 1
    end

    return n
end