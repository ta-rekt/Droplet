############# Interaction potential functions #############

# exports
export intPot!, getIntPot

# sets the self interaction potential of a species
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


############# Subspecies Interaction Functions #############
# MARK: -Subspecies

# These functions mirror the AbstractBulk pattern: all interacting Subspecies
# share a synchronized intSpecies list and intMat matrix.

# helper: find subspecies in intSpecies list, return its canonical index or nothing
function _findSubspeciesIdx(s::Subspecies, target::Subspecies)
    if s.intSpecies === nothing
        return nothing
    end
    for (sub, idx) in s.intSpecies
        if sub === target
            return idx
        end
    end
    return nothing
end

# helper: initialize interaction structures for a subspecies if not present
function _ensureIntStructures!(s::Subspecies)
    if s.intSpecies === nothing
        s.intSpecies = Vector{Tuple{Subspecies, Int}}()
    end
    if s.intMat === nothing
        s.intMat = Matrix{Float64}(undef, 0, 0)
    end
    if s.intIdx === nothing
        # subspecies not yet in its own list, add it
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

# helper: synchronize interaction matrix between two subspecies
# (they share the same matrix)
function _syncMatrices!(a::Subspecies, b::Subspecies)
    # both subspecies should point to the same matrix
    # merge their intSpecies lists and rebuild a shared matrix
    
    # collect all subspecies from both - id -> (subspecies, new_idx)
    allSubspecies = Dict{Int32, Tuple{Subspecies, Int}}()
    
    # first pass: collect unique subspecies and assign new indices
    newIdx = 1
    for (sub, _) in a.intSpecies
        if !haskey(allSubspecies, sub.id)
            allSubspecies[sub.id] = (sub, newIdx)
            newIdx += 1
        end
    end
    for (sub, _) in b.intSpecies
        if !haskey(allSubspecies, sub.id)
            allSubspecies[sub.id] = (sub, newIdx)
            newIdx += 1
        end
    end
    
    n = length(allSubspecies)
    newMat = zeros(Float64, n, n)
    
    # copy existing interaction values using old indices
    # from a's matrix
    for (sub1, oldIdx1) in a.intSpecies
        for (sub2, oldIdx2) in a.intSpecies
            (_, newIdx1) = allSubspecies[sub1.id]
            (_, newIdx2) = allSubspecies[sub2.id]
            if oldIdx1 <= size(a.intMat, 1) && oldIdx2 <= size(a.intMat, 2)
                newMat[newIdx1, newIdx2] = a.intMat[oldIdx1, oldIdx2]
            end
        end
    end
    
    # from b's matrix (only for subspecies not already covered)
    for (sub1, oldIdx1) in b.intSpecies
        for (sub2, oldIdx2) in b.intSpecies
            (_, newIdx1) = allSubspecies[sub1.id]
            (_, newIdx2) = allSubspecies[sub2.id]
            # only copy if not already set (prioritize existing values)
            if (newMat[newIdx1, newIdx2] == 0.0 
                && oldIdx1 <= size(b.intMat, 1) 
                && oldIdx2 <= size(b.intMat, 2))
                newMat[newIdx1, newIdx2] = b.intMat[oldIdx1, oldIdx2]
            end
        end
    end
    
    # build new intSpecies list
    newIntSpecies = Vector{Tuple{Subspecies, Int}}(undef, n)
    for (sub, idx) in values(allSubspecies)
        newIntSpecies[idx] = (sub, idx)
    end
    
    # update all subspecies to point to the same matrix and list
    for (sub, idx) in newIntSpecies
        sub.intSpecies = newIntSpecies
        sub.intMat = newMat
        sub.intIdx = idx
    end
end

"""
    intPot!(s::Subspecies, pot::Float64)

Set the self-interaction potential for a subspecies.
"""
function intPot!(s::Subspecies, pot::Float64)
    _ensureIntStructures!(s)
    
    s.intPot = pot
    
    # update the interaction matrix's diagonal entry
    s.intMat[s.intIdx, s.intIdx] = pot
    
    return s
end

"""
    intPot!(a::Subspecies, b::Subspecies, pot::Float64)

Set the interaction potential between two subspecies.
Updates the interaction matrix and interacting subspecies list.
This mirrors the AbstractBulk pattern.
"""
function intPot!(a::Subspecies, b::Subspecies, pot::Float64)
    # ensure both subspecies have interaction structures initialized
    _ensureIntStructures!(a)
    _ensureIntStructures!(b)
    
    # check if they already share the same matrix
    if a.intMat !== b.intMat
        # need to merge their interaction networks
        _syncMatrices!(a, b)
    end
    
    # now a and b share the same intMat and intSpecies
    idxA = a.intIdx
    idxB = b.intIdx
    
    # set the interaction potential (symmetric)
    a.intMat[idxA, idxB] = pot
    a.intMat[idxB, idxA] = pot
    
    b.intMat = a.intMat
    
    return nothing
end

"""
    getIntPot(a::Subspecies, b::Subspecies)

Get the interaction potential between two subspecies.
"""
function getIntPot(a::Subspecies, b::Subspecies)
    if a.intMat === nothing || a.intIdx === nothing
        return nothing
    end
    
    idxB = _findSubspeciesIdx(a, b)
    if idxB === nothing
        return nothing
    end
    
    return a.intMat[a.intIdx, idxB]
end


############# Polymer-Subspecies Convenience Functions #############
# MARK: -PolySubspecies

"""
    intPot!(p::Poly, subspeciesIdx::Int, pot::Float64)

Set the self-interaction potential for a specific subspecies of a polymer.
This is a convenience wrapper that calls intPot!(subspecies, pot).
"""
function intPot!(p::Poly, subspeciesIdx::Int, pot::Float64)
    if p.subspecies === nothing
        error("Subspecies not initialized. Call initSubspecies!(p, nSubspecies) first.")
    end
    
    if subspeciesIdx < 1 || subspeciesIdx > length(p.subspecies)
        error("Invalid subspecies index: $subspeciesIdx")
    end
    
    intPot!(p.subspecies[subspeciesIdx], pot)
    return p
end

"""
    intPot!(p::Poly, subIdx1::Int, subIdx2::Int, pot::Float64)

Set the interaction potential between two subspecies of the same polymer.
This is a convenience wrapper that calls intPot!(subspecies_a, subspecies_b, pot).
"""
function intPot!(p::Poly, subIdx1::Int, subIdx2::Int, pot::Float64)
    if p.subspecies === nothing
        error("Subspecies not initialized. Call initSubspecies!(p, nSubspecies) first.")
    end
    
    nSub = length(p.subspecies)
    if subIdx1 < 1 || subIdx1 > nSub || subIdx2 < 1 || subIdx2 > nSub
        error("Invalid subspecies indices: ($subIdx1, $subIdx2)")
    end
    
    intPot!(p.subspecies[subIdx1], p.subspecies[subIdx2], pot)
    
    return p
end

"""
    intPot!(p::Poly, subspeciesIdx::Int, b::AbstractBulk, pot::Float64)

Set the interaction potential between a polymer subspecies and a bulk species.
Stores in p.bulkIntPot dictionary using (subspecies.intIdx, bulk.id) as key.
"""
function intPot!(p::Poly, subspeciesIdx::Int, b::AbstractBulk, pot::Float64)
    if p.subspecies === nothing
        error("Subspecies not initialized. Call initSubspecies!(p, nSubspecies) first.")
    end
    
    if subspeciesIdx < 1 || subspeciesIdx > length(p.subspecies)
        error("Invalid subspecies index: $subspeciesIdx")
    end
    
    sub = p.subspecies[subspeciesIdx]
    
    # Ensure subspecies has intIdx
    _ensureIntStructures!(sub)
    
    # Initialize bulkIntPot dict if needed
    if p.bulkIntPot === nothing
        p.bulkIntPot = Dict{Tuple{Int, Int32}, Float64}()
    end
    
    # Initialize bulkSpecies list if needed
    if p.bulkSpecies === nothing
        p.bulkSpecies = Vector{AbstractBulk}()
    end
    
    # Add bulk species to list if not present
    if !(b in p.bulkSpecies)
        push!(p.bulkSpecies, b)
    end
    
    # Store interaction: key = (subspecies.intIdx, bulk.id)
    p.bulkIntPot[(sub.intIdx, b.id)] = pot
    
    return p
end

"""
    intPot!(p::Poly, b::AbstractBulk, pot::Float64)

Set the interaction potential between ALL subspecies of a polymer and a bulk species.
"""
function intPot!(p::Poly, b::AbstractBulk, pot::Float64)
    if p.subspecies === nothing
        # No subspecies - just set partner
        p.partner = b
        return p
    end
    
    # Set interaction for all subspecies
    for i in 1:length(p.subspecies)
        intPot!(p, i, b, pot)
    end
    
    p.partner = b
    
    return p
end

"""
    getIntPot(p::Poly, subIdx1::Int, subIdx2::Int)

Get the interaction potential between two subspecies of the same polymer.
"""
function getIntPot(p::Poly, subIdx1::Int, subIdx2::Int)
    if p.subspecies === nothing
        return p.intPot  # fallback to polymer's default
    end
    
    nSub = length(p.subspecies)
    if subIdx1 < 1 || subIdx1 > nSub || subIdx2 < 1 || subIdx2 > nSub
        return p.intPot
    end
    
    pot = getIntPot(p.subspecies[subIdx1], p.subspecies[subIdx2])
    return pot !== nothing ? pot : p.intPot
end

"""
    getIntPot(p::Poly, subspeciesIdx::Int, b::AbstractBulk)

Get the interaction potential between a polymer subspecies and a bulk species.
"""
function getIntPot(p::Poly, subspeciesIdx::Int, b::AbstractBulk)
    if p.subspecies === nothing || p.bulkIntPot === nothing
        return nothing
    end
    
    if subspeciesIdx < 1 || subspeciesIdx > length(p.subspecies)
        return nothing
    end
    
    sub = p.subspecies[subspeciesIdx]
    if sub.intIdx === nothing
        return nothing
    end
    
    key = (sub.intIdx, b.id)
    if haskey(p.bulkIntPot, key)
        return p.bulkIntPot[key]
    end
    
    return nothing
end
