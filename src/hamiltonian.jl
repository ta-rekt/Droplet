############# Polymer Energy functions #############
# MARK: -Polymer

# Hamiltonian of a polymer
export ham!
function ham!(p::Polymer)

    p.ham = 0
    totalChemEnergy = 0.0
    totalIntEnergy = 0.0
    
    key = p.head           # initial bead
    pos = p.pos[key][range(1, 3)]  # initial position

    for i in 1:p.length
        # Get potentials for this specific bead
        chemPot, intPot = getMonomerPotentials(p, Int32(i))
        
        # Chemical potential contribution
        totalChemEnergy -= chemPot
        
        # Interaction energy: neighbors - self
        # Note: nn(pos) counts occupied nns, so total_nn double-counts each interaction
        # Original formula: -intPot*N + intPot*total_nn
        # Per-bead: -intPot + intPot*nn_i
        nnn = nn(p, pos)
        totalIntEnergy += intPot * nnn - intPot
        
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
export deltaHam!
function deltaHam!(p::Polymer)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    p.deltaH = 0

    # bond addition
    if p.moveType == 1
        # Two new beads will be added at newPosi and newPosj
        # For now, assume they have the same potentials (uniform polymer case)
        # In a heterogeneous polymer, we'd need to know which monomer types these will be
        # This is a limitation of the current move proposal system
        chemPot, intPot = getMonomerPotentials(p, p.moveIdx)
        
        p.deltaH = (- 2 * (chemPot + intPot) 
                + intPot * (2*(nn(p, p.newPosi) + nn(p, p.newPosj)) + 2))

    # kink
    elseif p.moveType == 2

        if p.newPosi == p.bprev
            # bead ki is moving from bi to newPosj
            chemPot, intPot = getMonomerPotentials(p, p.moveIdx)
            p.deltaH = 2 * intPot * (nn(p, p.newPosj) - nn(p, p.bi))
        elseif p.newPosj == p.bnext
            # bead kj is moving from bj to newPosi
            chemPot, intPot = getMonomerPotentials(p, p.moveIdx + Int32(1))
            p.deltaH = 2 * intPot * (nn(p, p.newPosi) - nn(p, p.bj))
        else
            error("Move is not a kink")
        end

    # bond removal  
    elseif p.moveType == 3
        # Get potentials for the two beads being removed
        chemPoti, intPoti = getMonomerPotentials(p, p.moveIdx)
        chemPotj, intPotj = getMonomerPotentials(p, p.moveIdx + Int32(1))
        
        # When removing two bonded beads, the energy change is:
        # + chemical potential (removing negative contribution)
        # + 2*interaction potential (removing self-interaction, adjusted for the 
        # bond between them)
        # - interaction potential * neighbors
        p.deltaH = ((chemPoti + chemPotj) + 2*(intPoti + intPotj) 
                    - 2*intPoti * nn(p, p.bi) - 2*intPotj * nn(p, p.bj))
    end

    p.ham += p.deltaH

    return p
end


############# PolyMonomer energy functions #############
# MARK: -PolyMonomer

# Hamiltonian contribution from beads of this monomer type in a multicomponent polymer
# Note: PolyMonomer doesn't have its own ham field - energy is tracked at the Poly level
# This function calculates the contribution from beads of this specific monomer type
export ham!
function ham!(m::PolyMonomer)
    
    if m.intSpecies === nothing || m.intMat === nothing
        return 0.0
    end
    
    p = m.poly
    energy = 0.0
    
    # Count how many beads of this monomer type exist
    beadCount = 0
    
    # Iterate through polymer beads
    key = p.head
    beadIdx = 1
    
    for i in 1:p.length
        # Check if this bead belongs to this monomer type
        if p.monomerIndex[beadIdx] == findfirst(mon -> mon === m, p.monomers)
            beadCount += 1
            
            # Get position of this bead
            pos = p.pos[key][range(1, 3)]
            
            # Add chemical potential contribution
            energy -= m.chemPot
            
            # Add interaction contributions with all species
            for (sp, idx) in m.intSpecies
                selfFactor = m.intIdx == idx ? 0.5 : 1.0
                energy += m.intMat[m.intIdx, idx] * nnp(m, sp, pos) * selfFactor
            end
        end
        
        # Move to next bead
        key = p.pos[key][5]
        beadIdx += 1
        
        if key == 0
            break
        end
    end
    
    return energy
end


############# BulkMonomer energy functions #############
# MARK: -BulkMonomer

# Hamiltonian of a dilute species in the bulkMonomer
export ham!
function ham!(b::BulkMonomer)
    
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
export deltaHam!
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

# Hamiltonian of a dilute species in the bulkMonomer
export ham!
function ham!(b::FastBulkMonomer)
    
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
export deltaHam!
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
