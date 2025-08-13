# Hamiltonian of a polymer
export hamPoly!
function hamPoly!(s::Species)

    nnn = 0
    key = s.head           # initial bead
    pos = s.pos[key][1:3]  # initial posiion

    for i in 1:s.length
        nnn += nn(s, pos)
        key = s.pos[key][5]    # next bead
        pos = s.pos[key][1:3]  # next position
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    return - (s.chemPot + s.intPot) * s.length + s.intPot * nnn
end

# Number of nearest neighbors at a given position
export nn
function nn(s::Species, pos::Vector{Int32})
    n = 0

    for i in 1:3
        pos[i] += 1
        n += s.latS[pos...]

        pos[i] -= 2
        n += s.latS[pos...]

        pos[i] += 1
    end

    return n - 1
end