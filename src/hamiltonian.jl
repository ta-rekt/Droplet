
export ham!, deltaHam!

############# Polymer Energy functions #############
# MARK: -Polymer

# Hamiltonian of a polymer
function ham!(p::Polymer)
    p.ham = 0
    totalChemEnergy = 0.0
    totalIntEnergy = 0.0

    key = p.head           # initial bead
    pos = p.pos[key][range(1, 3)]  # initial position
    
     # assumes all monomers identical
    m = p.monomers[1]

    @assert !isnothing(m.intMat) "No interaction matrix!"

    for i in 1:p.length
        # Chemical potential contribution
        totalChemEnergy -= m.chemPot
        
        for (sp, idx) in m.intSpecies    
            # Interaction energy: neighbors - self, subtract 1 from self-interaction
            nnp_ = nnp(m, sp, pos)
            totalIntEnergy += m.intMat[m.intIdx, idx] * (nnp_ - Int(m.intIdx == idx))   
        end

        key = p.pos[key][5]    # next bead

        if key != 0
            pos = p.pos[key][range(1, 3)]  # next position
        end
    end

    @assert key == 0 begin
        "Tail is not N'th bead"
    end

    p.ham = totalChemEnergy + totalIntEnergy

    return p
end

# Energy increment of a proposed move. Requires proposeMove! to be called
# which sets moveType and newPosi/j. 
function deltaHam!(p::Polymer)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end
    
    p.deltaH = 0
    m = p.monomers[1]  # assumes all monomers identical
    
    @assert !isnothing(m.intMat) "No interaction matrix!"

    # bond addition
    if p.moveType == 1
        for (sp, idx) in m.intSpecies
            
            # bulk species interactions
            if m.intIdx != idx  
                sp.deltaHP = - 2 * (m.chemPot + m.intPot)  # change in chain length

                # bulk species energy increment
                sp.deltaHP += (m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosi) 
                             + m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosj))

                # polymer energy increment
                p.deltaH += sp.deltaHP
            
            # self-interactions
            else
                p.deltaH = - 2 * (m.chemPot + m.intPot)  # change in chain length

                p.deltaH += (2 * m.intPot * nn(p, p.newPosi) 
                           + 2 * m.intPot * nn(p, p.newPosj)
                           + 2 * m.intPot)  # beads i and j each gain a new neighbor
            end
        end

    # kink
    elseif p.moveType == 2
        for (sp, idx) in m.intSpecies

            # bead ki is moving from bi to newPosj
            if p.newPosi == p.bprev

                # bulk species interactions
                if m.intIdx != idx  

                    # bulk species energy increment
                    sp.deltaHP = (m.intMat[m.intIdx, idx] 
                                * (nnp(p, sp, p.newPosj) - nnp(p, sp, p.bi)))

                    # polymer energy increment
                    p.deltaH += sp.deltaHP

                # self-interactions
                else
                    p.deltaH = (2 * m.intMat[m.intIdx, idx] 
                                * (nnp(p, sp, p.newPosj) - nnp(p, sp, p.bi)))
                end

            # bead kj is moving from bj to newPosi
            elseif p.newPosj == p.bnext

                # bulk species interactions
                if m.intIdx != idx  

                    # bulk species energy increment
                    sp.deltaHP = (m.intMat[m.intIdx, idx] 
                                * (nnp(p, sp, p.newPosi) - nnp(p, sp, p.bj)))

                    # polymer energy increment
                    p.deltaH += sp.deltaHP

                # self-interactions
                else
                    p.deltaH = (2 * m.intMat[m.intIdx, idx] 
                                * (nnp(p, sp, p.newPosi) - nnp(p, sp, p.bj)))

                end

            else
                error("Move is not a kink")
            end
        end

    # bond removal  
    elseif p.moveType == 3
        for (sp, idx) in m.intSpecies
            
            # bulk species interactions
            if m.intIdx != idx  

                sp.deltaHP = 2 * (m.chemPot + m.intPot)   # change in chain length

                # bulk species energy increment
                sp.deltaHP -= (m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosi) 
                             + m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosj))

                # polymer energy increment
                p.deltaH += sp.deltaHP
            
            # self-interactions
            else
                p.deltaH = 2 * (m.chemPot + m.intPot)   # change in chain length

                p.deltaH -= (2 * m.intPot * nn(p, p.newPosi) 
                           + 2 * m.intPot * nn(p, p.newPosj)
                           - 2 * m.intPot) # beads i and j each lose a new neighbor
            end
        end
    end

    p.ham += p.deltaH

    return p
end

ham!(m::PolyMonomer) = ham!(m.poly)
deltaHam!(m::PolyMonomer) = deltaHam!(m.poly)

############# BulkMonomer energy functions #############
# MARK: -BulkMonomer

# Hamiltonian of a dilute species in the bulkMonomer
function ham!(b::AbstractBulkMonomer)
    
    # chemical potential
    b.ham = - b.chemPot * b.length

    # interaction potentials (including self)
    for (sp, idx) in b.intSpecies
        nnp_ = 0.0   # number of nearest partners
        selfFactor = b.intIdx == idx ? 0.5 : 1.0

        for i in 1:b.length
            pos = b.pos[i][range(1, 3)]
            nnp_ += nnp(b, sp, pos)
        end
        
        b.ham += b.intMat[b.intIdx, idx] * nnp_ * selfFactor
    end

    return b
end

# Energy increment for the bulkMonomer system
function deltaHam!(b::BulkMonomer)

    # energy increment from self-interactions
    b.deltaH = 0.0

    # energy increment from cross-interactions with each partner
    if !isnothing(b.intMat)
        for (sp, idx) in b.intSpecies
            sp.deltaHP = (b.intMat[b.intIdx, idx] 
                       * ((nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))
                       - Int(b.intIdx == idx))) # subtract 1 in self-interaction
            b.deltaH += sp.deltaHP
        end
    end

    b.ham += b.deltaH

    return b
end


############# FastBulkMonomer energy functions #############
# MARK: -FastBulkMonomer

# Energy increment for the bulkMonomer system
function deltaHam!(b::FastBulkMonomer)

    # energy increment from self-interactions
    b.deltaH = 0.0

    if !xor(b.moveIdx == 0, b.moveIdxNew == 0)
        return b
    else
        sign = Int(b.moveIdxNew == 0) * 2 - 1
    end

    # energy increment from cross-interactions with each partner
    if !isnothing(b.intMat)
        for (sp, idx) in b.intSpecies
            sp.deltaHP = (b.intMat[b.intIdx, idx] 
                       * (nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))) * sign
            b.deltaH += sp.deltaHP
        end
    end

    b.ham += b.deltaH

    return b
end
