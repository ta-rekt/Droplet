
export ham!, deltaHam!

############# Polymer Energy functions #############
# MARK: -Polymer

# Hamiltonian of a polymer
function ham!(p::Polymer)
    p.ham = 0
    totalChemEnergy = 0.0
    totalIntEnergy = 0.0

    key = p.head           # initial bead

    # Use first monomer to check intMat exists (all monomers share the same intMat)
    @assert !isnothing(p.monomers[1].intMat) "No interaction matrix!"

    for i in 1:p.length
        pos = p.pos[key][range(1, 3)]  # position of this bead
        m = p.monomers[p.pos[key][6]]  # get monomer type from pos[key][6]
        
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

        key = p.pos[key][5]    # next bead
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
    p.deltaHP = 0
    m = p.monomers[1]  # assumes all monomers identical
    
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
        
        # Energy change from type swaps at pattern tail (type X -> type 1)
        if currentType >= 2 && p.patternTails !== nothing && haskey(p.patternTails, currentType)
            _addSwapEnergyChange!(p, currentType, Int32(1), true)
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
        
        # Energy change from type swaps after pattern tail (type 1 -> type X)
        if currentType >= 2 && p.patternTails !== nothing && haskey(p.patternTails, currentType)
            _addSwapEnergyChange!(p, Int32(1), currentType, false)
        end
    end

    p.ham += p.deltaH

    return p
end

"""
    _addSwapEnergyChange!(p::Polymer, fromType::Int32, toType::Int32, isAddition::Bool)

Add energy change from type swaps to p.deltaH.
For addition: swaps happen at pattern tail (last 2 beads of pattern type change to base type)
For removal: swaps happen after pattern tail (next 2 base-type beads change to pattern type)
"""
function _addSwapEnergyChange!(p::Polymer, fromType::Int32, toType::Int32, isAddition::Bool)
    currentType = isAddition ? fromType : toType  # the pattern type
    tailKey = p.patternTails[currentType]
    
    # Get the two beads that will be swapped
    if isAddition
        # For addition: swap the tail and the bead before it
        swapKey1 = tailKey
        swapKey2 = p.pos[tailKey][4]
    else
        # For removal: swap the bead after tail and the one after that
        swapKey1 = p.pos[tailKey][5]
        swapKey2 = p.pos[swapKey1][5]
    end
    
    swapPos1 = p.pos[swapKey1][1:3]
    swapPos2 = p.pos[swapKey2][1:3]
    
    mOld = p.monomers[fromType]
    mNew = p.monomers[toType]
    
    # Calculate energy change for each swapped bead
    for (sp, idx) in mOld.intSpecies
        # Energy change = new interactions - old interactions
        # For each bead, subtract old type's interaction and add new type's interaction
        
        nnp1 = nnp(mOld, sp, swapPos1)
        nnp2 = nnp(mOld, sp, swapPos2)
        
        # Self-interactions need special handling
        if mOld.intIdx == idx
            # Old self-interactions being removed (with factor of 1/2 for each)
            p.deltaH -= mOld.intMat[mOld.intIdx, idx] * (nnp1 + nnp2) / 2.0
        else
            # Old cross-interactions being removed
            p.deltaH -= mOld.intMat[mOld.intIdx, idx] * (nnp1 + nnp2)
        end
    end
    
    for (sp, idx) in mNew.intSpecies
        nnp1 = nnp(mNew, sp, swapPos1)
        nnp2 = nnp(mNew, sp, swapPos2)
        
        if mNew.intIdx == idx
            # New self-interactions being added
            p.deltaH += mNew.intMat[mNew.intIdx, idx] * (nnp1 + nnp2) / 2.0
        else
            # New cross-interactions being added
            p.deltaH += mNew.intMat[mNew.intIdx, idx] * (nnp1 + nnp2)
        end
    end
    
    # Chemical potential change: 2 beads change type
    p.deltaH += 2 * mOld.chemPot - 2 * mNew.chemPot
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
