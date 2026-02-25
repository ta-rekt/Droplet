
export intPot!, intPot!, getIntPot, couplingPot!, getCouplingPot

############# Interaction potential functions #############

############# Helpers #############
# MARK: -Helpers

# helper: find species in intSpecies list, return its canonical index or nothing
function _findSpeciesIdx(s::AbstractMonomer, target::AbstractMonomer)
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
function _ensureIntStructures!(s::AbstractMonomer)
    if s.intSpecies === nothing
        s.intSpecies = Vector{Tuple{AbstractMonomer, Int}}()
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
function _syncMatrices!(a::AbstractMonomer, b::AbstractMonomer)
    # both species should point to the same matrix
    # merge their intSpecies lists and rebuild a shared matrix
    
    # collect all species from both - id -> (species, new_idx)
    allSpecies = Dict{Int, Tuple{AbstractMonomer, Int}}()
    
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
    newIntSpecies = Vector{Tuple{AbstractMonomer, Int}}(undef, n)
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


############# Monomers #############
# MARK: -Monomer

# sets the self interaction potential of a species
function intPot!(b::AbstractMonomer, pot::Float64)
    _ensureIntStructures!(b)
    
    b.intPot = pot

    # update the interaction matrix's diagonal entry
    b.intMat[b.intIdx, b.intIdx] = pot

    return b
end

# sets the interaction potential between two species
# updates the interaction matrix and interacting species list
function intPot!(a::AbstractMonomer, b::AbstractMonomer, pot::Float64)
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
function getIntPot(a::AbstractMonomer, b::AbstractMonomer)
    if a.intMat === nothing || a.intIdx === nothing
        return nothing
    end
    
    idxB = _findSpeciesIdx(a, b)
    if idxB === nothing
        return nothing
    end
    
    return a.intMat[a.intIdx, idxB]
end


############# Three-way coupling functions #############
# MARK: -Coupling

# helper: initialize coupling structures for a BulkMonomer if not present
function _ensureCouplingStructures!(s::BulkMonomer)
    # first ensure two-way interaction structures exist (needed for intIdx)
    _ensureIntStructures!(s)
    
    if s.couplingSpecies === nothing
        s.couplingSpecies = Vector{Tuple{AbstractBulkMonomer, Int32}}()
    end
    if s.couplingIdx === nothing
        # add self to coupling species list
        m = length(s.couplingSpecies) + 1
        s.couplingIdx = Int32(m)
        push!(s.couplingSpecies, (s, Int32(m)))
    end
    if s.couplingMat === nothing
        # initialize with dimensions: m (BulkMonomers) × n (all species) × n
        m = length(s.couplingSpecies)
        n = s.intMat !== nothing ? size(s.intMat, 1) : 1
        s.couplingMat = zeros(Float64, m, n, n)
    end
end

# helper: synchronize coupling structures between two BulkMonomers
function _syncCouplingStructures!(a::BulkMonomer, b::BulkMonomer)
    # both must already have coupling structures initialized
    # merge their couplingSpecies lists and rebuild shared couplingMat
    
    # first ensure intMat is synced (needed for partner/context indices)
    if a.intMat !== b.intMat
        _syncMatrices!(a, b)
    end
    
    # collect all BulkMonomers from both coupling lists
    allCouplingSpecies = Dict{Int32, Tuple{AbstractBulkMonomer, Int32}}()
    
    newIdx = Int32(1)
    for (sp, _) in a.couplingSpecies
        if !haskey(allCouplingSpecies, sp.id)
            allCouplingSpecies[sp.id] = (sp, newIdx)
            newIdx += Int32(1)
        end
    end
    for (sp, _) in b.couplingSpecies
        if !haskey(allCouplingSpecies, sp.id)
            allCouplingSpecies[sp.id] = (sp, newIdx)
            newIdx += Int32(1)
        end
    end
    
    m = length(allCouplingSpecies)
    n = size(a.intMat, 1)  # number of all species (from synced intMat)
    newCouplingMat = zeros(Float64, m, n, n)
    
    # copy existing coupling values from a's matrix
    for (sp, oldCoupIdx) in a.couplingSpecies
        (_, newCoupIdx) = allCouplingSpecies[sp.id]
        if oldCoupIdx <= size(a.couplingMat, 1)
            # copy all j,k values for this self species
            jmax = min(size(a.couplingMat, 2), n)
            kmax = min(size(a.couplingMat, 3), n)
            newCouplingMat[newCoupIdx, 1:jmax, 1:kmax] .= a.couplingMat[oldCoupIdx, 1:jmax, 1:kmax]
        end
    end
    
    # copy existing coupling values from b's matrix (if not already set)
    for (sp, oldCoupIdx) in b.couplingSpecies
        (_, newCoupIdx) = allCouplingSpecies[sp.id]
        if oldCoupIdx <= size(b.couplingMat, 1)
            jmax = min(size(b.couplingMat, 2), n)
            kmax = min(size(b.couplingMat, 3), n)
            for j in 1:jmax
                for k in 1:kmax
                    if newCouplingMat[newCoupIdx, j, k] == 0.0
                        newCouplingMat[newCoupIdx, j, k] = b.couplingMat[oldCoupIdx, j, k]
                    end
                end
            end
        end
    end
    
    # build new couplingSpecies list
    newCouplingSpecies = Vector{Tuple{AbstractBulkMonomer, Int32}}(undef, m)
    for (sp, idx) in values(allCouplingSpecies)
        newCouplingSpecies[idx] = (sp, idx)
    end
    
    # update all BulkMonomers to point to same structures
    for (sp, idx) in newCouplingSpecies
        sp.couplingSpecies = newCouplingSpecies
        sp.couplingMat = newCouplingMat
        sp.couplingIdx = idx
    end
end

# helper: expand coupling matrix when intMat grows (new species added)
function _expandCouplingMat!(s::BulkMonomer)
    if s.couplingMat === nothing
        return
    end
    
    n_old = size(s.couplingMat, 2)
    n_new = size(s.intMat, 1)
    
    if n_new > n_old
        m = size(s.couplingMat, 1)
        newCouplingMat = zeros(Float64, m, n_new, n_new)
        newCouplingMat[:, 1:n_old, 1:n_old] .= s.couplingMat
        s.couplingMat = newCouplingMat
        
        # update all BulkMonomers sharing this coupling matrix
        if s.couplingSpecies !== nothing
            for (sp, _) in s.couplingSpecies
                sp.couplingMat = newCouplingMat
            end
        end
    end
end

"""
    couplingPot!(self::BulkMonomer, partner::AbstractMonomer, context::AbstractMonomer, pot::Float64)

Set the three-way coupling potential for a BulkMonomer.

The coupling contributes to `self`'s Hamiltonian when:
- `self` occupies the same lattice site as `context`
- `self` has a nearest neighbor of type `partner`

The coupling is symmetric in partner/context: setting (self, A, B, pot) also sets (self, B, A, pot).

All three species (self, partner, context) must be distinct.

# Arguments
- `self::BulkMonomer`: The BulkMonomer whose Hamiltonian is modified (cannot be FastBulkMonomer)
- `partner::AbstractMonomer`: The species that must be a nearest neighbor
- `context::AbstractMonomer`: The species that must occupy the same site as self
- `pot::Float64`: The coupling strength (energy contribution per neighbor)

# Example
```julia
a = BulkMonomer(; name="A", length=100)
b = BulkMonomer(; name="B", length=100)
c = BulkMonomer(; name="C", length=100)

# Set up standard two-way interactions first
intPot!(a, -1.0)
intPot!(b, -1.0)
intPot!(c, -1.0)

# When A is at the same site as C, A's interaction with B neighbors is modified by -0.5
couplingPot!(a, b, c, -0.5)
```
"""
function couplingPot!(self::BulkMonomer, partner::AbstractMonomer, context::AbstractMonomer, pot::Float64)
    # validate: all three must be distinct
    if self === partner || self === context || partner === context
        error("couplingPot!: self, partner, and context must all be distinct species")
    end
    
    # ensure two-way interaction structures exist for all species
    _ensureIntStructures!(self)
    _ensureIntStructures!(partner)
    _ensureIntStructures!(context)
    
    # sync two-way matrices so all species share the same intMat
    if self.intMat !== partner.intMat
        _syncMatrices!(self, partner)
    end
    if self.intMat !== context.intMat
        _syncMatrices!(self, context)
    end
    
    # initialize coupling structures for self
    _ensureCouplingStructures!(self)
    
    # expand coupling matrix if intMat has grown
    _expandCouplingMat!(self)
    
    # get indices
    selfCoupIdx = self.couplingIdx
    partnerIntIdx = partner.intIdx
    contextIntIdx = context.intIdx
    
    # set coupling (symmetric in partner/context)
    self.couplingMat[selfCoupIdx, partnerIntIdx, contextIntIdx] = pot
    self.couplingMat[selfCoupIdx, contextIntIdx, partnerIntIdx] = pot
    
    return nothing
end

# get coupling potential for a three-way interaction
function getCouplingPot(self::BulkMonomer, partner::AbstractMonomer, context::AbstractMonomer)
    if self.couplingMat === nothing || self.couplingIdx === nothing
        return nothing
    end
    if partner.intIdx === nothing || context.intIdx === nothing
        return nothing
    end
    
    # check bounds
    if partner.intIdx > size(self.couplingMat, 2) || context.intIdx > size(self.couplingMat, 3)
        return nothing
    end
    
    return self.couplingMat[self.couplingIdx, partner.intIdx, context.intIdx]
end