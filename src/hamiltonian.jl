############# Polymer Energy functions #############

# Hamiltonian of a polymer
export ham!
function ham!(p::Poly)

    p.ham = 0
    nnn = 0                # number of nearest neighbors
    key = p.head           # initial bead
    pos = p.pos[key][range(1, 3)]  # initial posiion

    for i in 1:p.length
        nnn += nn(p, pos)
        key = p.pos[key][5]    # next bead

        if key != 0
            pos = p.pos[key][range(1, 3)]  # next position
        end
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    p.ham = - (p.chemPot + p.intPot) * p.length + p.intPot * nnn

    return p
end

# Energy increment of a proposed move. Requires proposeMove! to be called
# which sets moveType and newPosi/j. 
export deltaHam!
function deltaHam!(p::Poly)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    p.deltaH = 0

    # bond addition
    if p.moveType == 1
        p.deltaH = (- 2 * (p.chemPot + p.intPot) 
                    + p.intPot * (2*(nn(p, p.newPosi) + nn(p, p.newPosj)) + 2))

    # kink
    elseif p.moveType == 2

        if p.newPosi == p.bprev 
            p.deltaH = 2 * p.intPot * (nn(p, p.newPosj) - nn(p, p.bi))
        elseif p.newPosj == p.bnext
            p.deltaH = 2 * p.intPot * (nn(p, p.newPosi) - nn(p, p.bj))
        else
            error("Move is not a kink")
        end

    # bond removal  
    elseif p.moveType == 3 
        p.deltaH = (2 * (p.chemPot + p.intPot)  
                    - p.intPot * (2*(nn(p, p.bi) + nn(p, p.bj)) - 2))
    end

    p.ham += p.deltaH

    return p
end


############# Bulk energy functions #############

# Hamiltonian of a dilute species in the bulk
export ham!
function ham!(b::AbstractBulk)
    
    b.ham = 0
    nnn = 0    # number of nearest neighbors
    nnnp = 0   # number of nearest neighbors with the partner species

    for i in 1:b.length
        pos = b.pos[i][range(1, 3)]
        nnn += nn(b, pos)
    end

    if isnothing(b.partner)
        b.intPotP = 0
    else
        p = b.partner

        @assert !isnothing(b.intPotP) "No inter-partner interaction potential!"
        @assert !(b.intPotP == p.intPotP) "Non-reciprocal interaction!"

        for i in 1:b.length
            pos = b.pos[i][range(1, 3)]
            nnnp += nnp(b, pos)
        end
    end

    b.ham = - b.chemPot * b.length + b.intPot * nnn + b.intPotP * nnnp

    return b
end

# Energy increment for the bulk system
export deltaHam!
function deltaHam!(b::AbstractBulk)

    b.deltaH = b.intPot * (2*(nn(b, b.newPos) - nn(b, b.bi)) - 2)

    if !isnothing(b.partner)

        @assert !isnothing(b.intPotP) "No inter-partner interaction potential!"
        @assert !(b.intPotP == p.intPotP) "Non-reciprocal interaction!"

        b.deltaH += b.intPotP * (2*(nnp(b, b.newPos) - nnp(b, b.bi)) - 2)
    end

    b.ham += b.deltaH

    return b
end
