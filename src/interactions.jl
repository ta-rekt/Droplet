############# Interaction potential functions #############

# sets the self interaction potential of a species
export intPot!
function intPot!(b::AbstractBulk, pot::Float64)
    _ensureIntStructures!(b)
    
    b.intPot = pot

    # update the interaction matrix's diagonal entry
    b.intMat[b.intIdx, b.intIdx] = pot

    return b
end

# helper: find species in intSpecies list, return its canonical index or nothing
function _findSpeciesIdx(s::AbstractBulk, target::AbstractBulk)
    if s.intSpecies === nothing
        return nothing
    end
    for (sp, idx) in s.intSpecies
        if sp === target
            return idx
        end
    end
    return nothing
end

# helper: initialize interaction structures for a species if not present
function _ensureIntStructures!(s::AbstractBulk)
    if s.intSpecies === nothing
        s.intSpecies = Vector{Tuple{AbstractBulk, Int}}()
    end
    if s.intMat === nothing
        s.intMat = Matrix{Float64}(undef, 0, 0)
    end
    if s.intIdx === nothing
        # species not yet in its own list, add it
        n = length(s.intSpecies) + 1
        s.intIdx = n
        push!(s.intSpecies, (s, n))
        
        # expand matrix to include self
        newMat = zeros(Float64, n, n)
        if n > 1
            newMat[1:n-1, 1:n-1] .= s.intMat
        end
        newMat[n, n] = s.intPot !== nothing ? s.intPot : 0.0
        s.intMat = newMat
    end
end

# helper: synchronize interaction matrix between two species 
# (they share the same matrix)
function _syncMatrices!(a::AbstractBulk, b::AbstractBulk)
    # both species should point to the same matrix
    # merge their intSpecies lists and rebuild a shared matrix
    
    # collect all species from both - id -> (species, new_idx)
    allSpecies = Dict{Int, Tuple{AbstractBulk, Int}}() 
    
    # first pass: collect unique species and assign new indices
    newIdx = 1
    for (sp, _) in a.intSpecies
        if !haskey(allSpecies, sp.id)
            allSpecies[sp.id] = (sp, newIdx)
            newIdx += 1
        end
    end
    for (sp, _) in b.intSpecies
        if !haskey(allSpecies, sp.id)
            allSpecies[sp.id] = (sp, newIdx)
            newIdx += 1
        end
    end
    
    n = length(allSpecies)
    newMat = zeros(Float64, n, n)
    
    # copy existing interaction values using old indices
    # from a's matrix
    for (sp1, oldIdx1) in a.intSpecies
        for (sp2, oldIdx2) in a.intSpecies
            (_, newIdx1) = allSpecies[sp1.id]
            (_, newIdx2) = allSpecies[sp2.id]
            if oldIdx1 <= size(a.intMat, 1) && oldIdx2 <= size(a.intMat, 2)
                newMat[newIdx1, newIdx2] = a.intMat[oldIdx1, oldIdx2]
            end
        end
    end
    
    # from b's matrix (only for species not already covered)
    for (sp1, oldIdx1) in b.intSpecies
        for (sp2, oldIdx2) in b.intSpecies
            (_, newIdx1) = allSpecies[sp1.id]
            (_, newIdx2) = allSpecies[sp2.id]
            # only copy if not already set (prioritize existing values)
            if (newMat[newIdx1, newIdx2] == 0.0 
                && oldIdx1 <= size(b.intMat, 1) 
                && oldIdx2 <= size(b.intMat, 2))
                newMat[newIdx1, newIdx2] = b.intMat[oldIdx1, oldIdx2]
            end
        end
    end
    
    # build new intSpecies list
    newIntSpecies = Vector{Tuple{AbstractBulk, Int}}(undef, n)
    for (sp, idx) in values(allSpecies)
        newIntSpecies[idx] = (sp, idx)
    end
    
    # update all species to point to the same matrix and list
    for (sp, idx) in newIntSpecies
        sp.intSpecies = newIntSpecies
        sp.intMat = newMat
        sp.intIdx = idx
    end
end

# sets the interaction potential between two species
# updates the interaction matrix and interacting species list
export intPot!
function intPot!(a::AbstractBulk, b::AbstractBulk, pot::Float64)
    # ensure both species have interaction structures initialized
    _ensureIntStructures!(a)
    _ensureIntStructures!(b)
    
    # check if they already share the same matrix
    if a.intMat !== b.intMat
        # need to merge their interaction networks
        _syncMatrices!(a, b)
    end
    
    # now a and b share the same intMat and intSpecies
    # check if b is already in a's list (should be after sync)
    idxA = a.intIdx
    idxB = b.intIdx
    
    # set the interaction potential (symmetric)
    a.intMat[idxA, idxB] = pot
    a.intMat[idxB, idxA] = pot

    b.intMat = a.intMat
    
    return nothing
end

# get interaction potential between two species
export getIntPot
function getIntPot(a::AbstractBulk, b::AbstractBulk)
    if a.intMat === nothing || a.intIdx === nothing
        return nothing
    end
    
    idxB = _findSpeciesIdx(a, b)
    if idxB === nothing
        return nothing
    end
    
    return a.intMat[a.intIdx, idxB]
end
