
export ham!, deltaHam!

############# Polymer Energy functions #############
# MARK: -Polymer

# Hamiltonian of a polymer
function ham!(p::Polymer)
    p.ham = 0
    totalChemEnergy = 0.0
    totalIntEnergy = 0.0

    key = p.head
    
    # Get any monomer to check intMat exists (all share same matrix)
    refMonomer = first(values(p.monomers))
    @assert !isnothing(refMonomer.intMat) "No interaction matrix!"

    for i in 1:p.length
        pos = p.pos[key][1:3]
        m = p.monomers[p.pos[key][6]]  # Dictionary lookup by intIdx
        
        # Chemical potential contribution
        totalChemEnergy -= m.chemPot
        
        for (sp, idx) in m.intSpecies    
            # Interaction energy: count neighbors
            nnp_ = nnp(m, sp, pos)
            
            # For self-interactions, divide by 2 to avoid double-counting bonds
            # For cross-interactions, count all neighbors
            if m.intIdx == idx
                totalIntEnergy += m.intMat[m.intIdx, idx] * nnp_ / 2.0
            else
                totalIntEnergy += m.intMat[m.intIdx, idx] * nnp_
            end
        end

        key = p.pos[key][5]
    end

    @assert key == 0 "Tail is not N'th bead"

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
    p.deltaHP = 0
    
    # Get monomer for the bead being moved (ki)
    m = p.monomers[p.pos[p.ki][6]]
    
    @assert !isnothing(m.intMat) "No interaction matrix!"

    # bond addition
    if p.moveType == 1
        for (sp, idx) in m.intSpecies
            
            # bulk species interactions
            if m.intIdx != idx  

                # bulk species energy increment
                setDeltaHP(sp, (m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosi) 
                              + m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosj)))

                # polymer energy increment
                p.deltaH += getDeltaHP(sp)
            
            # self-interactions
            else
                # Chemical potential: adding 2 beads
                p.deltaH += - 2 * m.chemPot
                
                # Interaction energy: count bonds being added
                # New beads at newPosi and newPosj will be connected (1 bond)
                # New bead at newPosi will connect to nn(newPosi) existing neighbors
                # New bead at newPosj will connect to nn(newPosj) existing neighbors  
                # Total bonds added: 1 + nn(newPosi) + nn(newPosj)
                p.deltaH += m.intPot * (nn(p, p.newPosi) + nn(p, p.newPosj) + 1)
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
                    setDeltaHP(sp, (m.intMat[m.intIdx, idx] 
                                 * (nnp(p, sp, p.newPosj) - nnp(p, sp, p.bi))))

                    # polymer energy increment
                    p.deltaH += getDeltaHP(sp)

                # self-interactions
                else
                    # Bead ki moves from bi to newPosj
                    # Energy change = (interactions at newPosj) - (interactions at bi)
                    # With division by 2: each interaction counted once per bond
                    p.deltaH += m.intPot * (nn(p, p.newPosj) - nn(p, p.bi))
                end

            # bead kj is moving from bj to newPosi
            elseif p.newPosj == p.bnext

                # bulk species interactions
                if m.intIdx != idx  

                    # bulk species energy increment
                    setDeltaHP(sp, (m.intMat[m.intIdx, idx] 
                                 * (nnp(p, sp, p.newPosi) - nnp(p, sp, p.bj))))

                    # polymer energy increment
                    p.deltaH += getDeltaHP(sp)

                # self-interactions
                else
                    # Bead kj moves from bj to newPosi
                    # Energy change = (interactions at newPosi) - (interactions at bj)
                    p.deltaH += m.intPot * (nn(p, p.newPosi) - nn(p, p.bj))
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

                # bulk species energy increment
                # Use bi and bj (positions being removed), not newPosi/newPosj
                setDeltaHP(sp, -(m.intMat[m.intIdx, idx] * nnp(p, sp, p.bi) 
                               + m.intMat[m.intIdx, idx] * nnp(p, sp, p.bj)))

                # polymer energy increment
                p.deltaH += getDeltaHP(sp)
            
            # self-interactions
            else
                # Chemical potential: removing 2 beads
                p.deltaH += 2 * m.chemPot
                
                # Interaction energy: count bonds being removed
                # Beads bi and bj are connected (1 bond between them)
                # Bead bi has nn(bi) - 1 other neighbors (excluding bj)
                # Bead bj has nn(bj) - 1 other neighbors (excluding bi)
                # Total bonds removed: 1 + (nn(bi)-1) + (nn(bj)-1) = nn(bi) + nn(bj) - 1
                p.deltaH -= m.intPot * (nn(p, p.bi) + nn(p, p.bj) - 1)
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
            setDeltaHP(sp, (b.intMat[b.intIdx, idx] 
                       * ((nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))
                       - Int(b.intIdx == idx)))) # subtract 1 in self-interaction
            b.deltaH += getDeltaHP(sp)
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
            setDeltaHP(sp, (b.intMat[b.intIdx, idx] 
                       * (nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))) * sign)
            b.deltaH += getDeltaHP(sp)
        end
    end

    b.ham += b.deltaH

    return b
end
