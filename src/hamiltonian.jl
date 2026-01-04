############# Polymer Energy functions #############
# MARK: -Poly

export ham!, deltaHam!

"""
    ham!(p::Poly)

Calculate the Hamiltonian of a polymer, accounting for bead-dependent interactions
via subspecies if defined.

When subspecies are defined:
- Each bead contributes its subspecies' chemPot
- Self-interactions use the subspecies' intPot when neighbors are same subspecies
- Cross-subspecies interactions use the defined cross-interaction potential

When subspecies are NOT defined (backward compatible):
- Uses p.chemPot and p.intPot uniformly for all beads
"""
function ham!(p::Poly)
    p.ham = 0.0
    
    # Check if using subspecies-based interactions
    if p.subspecies !== nothing && p.beadTypes !== nothing
        ham_subspecies!(p)
    else
        ham_legacy!(p)
    end
    
    return p
end

# Legacy hamiltonian (backward compatible)
function ham_legacy!(p::Poly)
    nnn = 0                # number of nearest neighbors
    key = p.head           # initial bead
    pos = p.pos[key][range(1, 3)]  # initial position

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

# Subspecies-aware hamiltonian
function ham_subspecies!(p::Poly)
    p.ham = 0.0
    
    key = p.head
    beadIdx = 1
    
    while key != 0
        subIdx = p.beadTypes[beadIdx]
        sub = p.subspecies[subIdx]
        pos = p.pos[key][range(1, 3)]
        
        # Chemical potential contribution
        p.ham -= sub.chemPot !== nothing ? sub.chemPot : p.chemPot
        
        # Self-interaction contribution (with same subspecies neighbors)
        # We count nn here for lattice neighbors of same subspecies
        intPot_self = sub.intPot !== nothing ? sub.intPot : p.intPot
        p.ham -= intPot_self  # bonding energy
        nnn_same = nn(p, pos)  # neighbors on lattice
        p.ham += intPot_self * nnn_same
        
        key = p.pos[key][5]
        beadIdx += 1
    end
    
    # Add bulk species interactions if defined
    if p.bulkSpecies !== nothing && p.bulkIntPot !== nothing
        for bulk in p.bulkSpecies
            key = p.head
            beadIdx = 1
            while key != 0
                subIdx = p.beadTypes[beadIdx]
                sub = p.subspecies[subIdx]
                pos = p.pos[key][range(1, 3)]
                
                # Look up interaction using subspecies' intIdx
                if sub.intIdx !== nothing
                    intKey = (sub.intIdx, bulk.id)
                    if haskey(p.bulkIntPot, intKey)
                        intPot_sp = p.bulkIntPot[intKey]
                        p.ham += intPot_sp * nnp(p, bulk, pos)
                    end
                end
                
                key = p.pos[key][5]
                beadIdx += 1
            end
        end
    elseif p.partner !== nothing
        # Legacy partner interaction (no subspecies-specific potential)
        key = p.head
        while key != 0
            pos = p.pos[key][range(1, 3)]
            # Use polymer's intPot as default for partner interaction
            p.ham += p.intPot * nnp(p, p.partner, pos)
            key = p.pos[key][5]
        end
    end
    
    return p
end


"""
    deltaHam!(p::Poly)

Energy increment of a proposed move, accounting for bead-dependent interactions.
Requires proposeMove! to be called which sets moveType and newPosi/j.
"""
function deltaHam!(p::Poly)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    p.deltaH = 0.0
    
    # Check if using subspecies-based interactions
    if p.subspecies !== nothing && p.beadTypes !== nothing
        deltaHam_subspecies!(p)
    else
        deltaHam_legacy!(p)
    end

    p.ham += p.deltaH

    return p
end

# Legacy deltaHam (backward compatible)
function deltaHam_legacy!(p::Poly)
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
    
    return p
end

# Subspecies-aware deltaHam
function deltaHam_subspecies!(p::Poly)
    # Get subspecies indices for the involved beads
    # ki and kj are the keys for the two beads involved in the move
    subIdx_i = getSubspeciesIdx(p, p.ki)
    subIdx_j = getSubspeciesIdx(p, p.kj)
    
    sub_i = p.subspecies[subIdx_i]
    sub_j = p.subspecies[subIdx_j]
    
    # Get interaction potentials
    chemPot_i = sub_i.chemPot !== nothing ? sub_i.chemPot : p.chemPot
    chemPot_j = sub_j.chemPot !== nothing ? sub_j.chemPot : p.chemPot
    intPot_i = sub_i.intPot !== nothing ? sub_i.intPot : p.intPot
    intPot_j = sub_j.intPot !== nothing ? sub_j.intPot : p.intPot
    
    # Cross-subspecies interaction (for the bond between i and j)
    # Use the matrix-based getIntPot which uses Subspecies intMat
    intPot_ij = getIntPot(sub_i, sub_j)
    if intPot_ij === nothing
        intPot_ij = p.intPot  # fallback
    end
    
    # bond addition: adding two new beads
    if p.moveType == 1
        # Energy cost: chemical potential for two new beads
        p.deltaH = - (chemPot_i + intPot_i) - (chemPot_j + intPot_j)
        
        # Energy gain from neighbors (each new bead gains nn interactions)
        p.deltaH += intPot_i * nn(p, p.newPosi) + intPot_j * nn(p, p.newPosj)
        
        # Bond energy between the two new beads
        p.deltaH += intPot_ij

    # kink: moving one bead to a new position
    elseif p.moveType == 2
        if p.newPosi == p.bprev 
            # Moving bead j
            p.deltaH = intPot_j * (nn(p, p.newPosj) - nn(p, p.bj))
            # Account for bond with bead i (unchanged) and next bead
            p.deltaH += intPot_j * (nn(p, p.newPosj) - nn(p, p.bj))
        elseif p.newPosj == p.bnext
            # Moving bead i
            p.deltaH = intPot_i * (nn(p, p.newPosi) - nn(p, p.bi))
            p.deltaH += intPot_i * (nn(p, p.newPosi) - nn(p, p.bi))
        else
            error("Move is not a kink")
        end

    # bond removal: removing two beads
    elseif p.moveType == 3 
        # Energy gain: chemical potential returned
        p.deltaH = (chemPot_i + intPot_i) + (chemPot_j + intPot_j)
        
        # Energy cost from losing neighbors
        p.deltaH -= intPot_i * nn(p, p.bi) + intPot_j * nn(p, p.bj)
        
        # Bond energy between the two beads being removed
        p.deltaH -= intPot_ij
    end
    
    # Add bulk species interaction energy changes
    if p.bulkSpecies !== nothing && p.bulkIntPot !== nothing
        for bulk in p.bulkSpecies
            # Get interaction potentials for each subspecies with this bulk
            intPot_i_sp = nothing
            intPot_j_sp = nothing
            
            if sub_i.intIdx !== nothing
                intKey_i = (sub_i.intIdx, bulk.id)
                if haskey(p.bulkIntPot, intKey_i)
                    intPot_i_sp = p.bulkIntPot[intKey_i]
                end
            end
            
            if sub_j.intIdx !== nothing
                intKey_j = (sub_j.intIdx, bulk.id)
                if haskey(p.bulkIntPot, intKey_j)
                    intPot_j_sp = p.bulkIntPot[intKey_j]
                end
            end
            
            if p.moveType == 1  # bond addition
                if intPot_i_sp !== nothing
                    p.deltaH += intPot_i_sp * nnp(p, bulk, p.newPosi)
                end
                if intPot_j_sp !== nothing
                    p.deltaH += intPot_j_sp * nnp(p, bulk, p.newPosj)
                end
            elseif p.moveType == 2  # kink
                if p.newPosi == p.bprev && intPot_j_sp !== nothing
                    p.deltaH += intPot_j_sp * (nnp(p, bulk, p.newPosj) - nnp(p, bulk, p.bj))
                elseif p.newPosj == p.bnext && intPot_i_sp !== nothing
                    p.deltaH += intPot_i_sp * (nnp(p, bulk, p.newPosi) - nnp(p, bulk, p.bi))
                end
            elseif p.moveType == 3  # bond removal
                if intPot_i_sp !== nothing
                    p.deltaH -= intPot_i_sp * nnp(p, bulk, p.bi)
                end
                if intPot_j_sp !== nothing
                    p.deltaH -= intPot_j_sp * nnp(p, bulk, p.bj)
                end
            end
        end
    elseif p.partner !== nothing
        # Legacy partner interaction (no subspecies-specific)
        if p.moveType == 1
            p.deltaH += p.intPot * (nnp(p, p.partner, p.newPosi) + nnp(p, p.partner, p.newPosj))
        elseif p.moveType == 2
            if p.newPosi == p.bprev 
                p.deltaH += p.intPot * (nnp(p, p.partner, p.newPosj) - nnp(p, p.partner, p.bj))
            elseif p.newPosj == p.bnext
                p.deltaH += p.intPot * (nnp(p, p.partner, p.newPosi) - nnp(p, p.partner, p.bi))
            end
        elseif p.moveType == 3
            p.deltaH -= p.intPot * (nnp(p, p.partner, p.bi) + nnp(p, p.partner, p.bj))
        end
    end
    
    return p
end


############# Bulk energy functions #############
# MARK: -Bulk

# Hamiltonian of a dilute species in the bulk
function ham!(b::Bulk)
    
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

# Energy increment for the bulk system
function deltaHam!(b::Bulk)

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


############# FastBulk energy functions #############
# MARK: -FastBulk

# Hamiltonian of a dilute species in the bulk
function ham!(b::FastBulk)
    
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


# Energy increment for the bulk system
function deltaHam!(b::FastBulk)

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
