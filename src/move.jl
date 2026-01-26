
export  overlap, addBond!, kinkBond!, removeBond!, identifyBond!, 
        makeMove!, move!, identifyBond!

############# Polymer move functions #############
# MARK: -Polymer

function overlap(p::Polymer)
    return ((p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] == 1 && 
             !coords_equal(p.bprev, p.newPosi)) ||
            (p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] == 1 && 
             !coords_equal(p.bnext, p.newPosj)))
end

function addBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi != p.bprev) && (p.newPosj != p.bnext) begin
        "Move is not a bond addition"
    end

    # set new lattice sites to 1
    p.latS[p.newPosi...] = 1
    p.latS[p.newPosj...] = 1

    # Find next available keys
    next_key1 = KEYCOUNT[] + 1
    next_key2 = KEYCOUNT[] + 2
    while haskey(p.pos, next_key1) || haskey(p.pos, next_key2)
        KEYCOUNT[] += 1
        next_key1 = KEYCOUNT[] + 1
        next_key2 = KEYCOUNT[] + 2
    end

    # Get monomer type from parent bead
    currentType = p.pos[p.ki][6]

    # add two new entries for the bond
    p.pos[next_key1] = Vector{Int32}(undef, 6)
    p.pos[next_key1][1] = p.newPosi[1]
    p.pos[next_key1][2] = p.newPosi[2]
    p.pos[next_key1][3] = p.newPosi[3]
    p.pos[next_key1][4] = p.ki
    p.pos[next_key1][5] = next_key2
    p.pos[next_key1][6] = currentType  # inherit monomer type

    p.pos[next_key2] = Vector{Int32}(undef, 6)
    p.pos[next_key2][1] = p.newPosj[1]
    p.pos[next_key2][2] = p.newPosj[2]
    p.pos[next_key2][3] = p.newPosj[3]
    p.pos[next_key2][4] = next_key1
    p.pos[next_key2][5] = p.kj
    p.pos[next_key2][6] = currentType  # inherit monomer type

    # re-route the old pointers
    p.pos[p.ki][5] = next_key1
    p.pos[p.kj][4] = next_key2

    # Update keys vector
    push!(p.keys_vec, next_key1)
    push!(p.keys_vec, next_key2)

    KEYCOUNT[] += 2

    # If adding inside a pattern (type >= 2), swap tail beads to base type
    if currentType >= 2 && p.patternTails !== nothing && haskey(p.patternTails, currentType)
        tailKey = p.patternTails[currentType]
        prevTailKey = p.pos[tailKey][4]
        
        # Swap last two beads of pattern to base type 1
        p.pos[tailKey][6] = Int32(1)
        p.pos[prevTailKey][6] = Int32(1)
        
        # Update patternTails: new tail is two beads back
        newTailKey = p.pos[prevTailKey][4]
        if newTailKey != 0 && p.pos[newTailKey][6] == currentType
            p.patternTails[currentType] = newTailKey
        else
            # Pattern exhausted (shouldn't happen with proper validation)
            delete!(p.patternTails, currentType)
        end
    end

    return p
end

function kinkBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) ⊻ (p.newPosj == p.bnext) begin
        "Move is not a kink"
    end

    # set angled site to 0, flat site to 1 (allocation-free: no splatting or array construction)
    if p.newPosi == p.bprev 
        p.latS[p.bi[1], p.bi[2], p.bi[3]] = 0
        p.latS[p.newPosj[1], p.newPosj[2], p.newPosj[3]] = 1

        # Preserve monomer type when creating new entry
        oldType = p.pos[p.ki][6]
        p.pos[p.ki] = Vector{Int32}(undef, 6)
        p.pos[p.ki][1] = p.newPosj[1]
        p.pos[p.ki][2] = p.newPosj[2]
        p.pos[p.ki][3] = p.newPosj[3]
        p.pos[p.ki][4] = p.kprev
        p.pos[p.ki][5] = p.kj
        p.pos[p.ki][6] = oldType

    elseif p.newPosj == p.bnext
        p.latS[p.newPosi[1], p.newPosi[2], p.newPosi[3]] = 1
        p.latS[p.bj[1], p.bj[2], p.bj[3]] = 0

        # Preserve monomer type when creating new entry
        oldType = p.pos[p.kj][6]
        p.pos[p.kj] = Vector{Int32}(undef, 6)
        p.pos[p.kj][1] = p.newPosi[1]
        p.pos[p.kj][2] = p.newPosi[2]
        p.pos[p.kj][3] = p.newPosi[3]
        p.pos[p.kj][4] = p.ki
        p.pos[p.kj][5] = p.knext
        p.pos[p.kj][6] = oldType
    end

    return p
end

function removeBond!(p::Polymer)

    # TODO: handle corner cases
    @assert (p.newPosi == p.bprev) && (p.newPosj == p.bnext) begin
        "Move is not a bond removal"
    end
    
    # Get current type before removal for pattern handling
    currentType = p.pos[p.ki][6]
    
    # reset old lattice sites to 0
    p.latS[p.bi...] = 0
    p.latS[p.bj...] = 0

    # join ends together - handle boundary cases
    if p.kprev != 0
        p.pos[p.kprev][5] = p.knext
    else
        # Removing from the start, update head
        p.head = p.knext
    end
    
    if p.knext != 0
        p.pos[p.knext][4] = p.kprev
    end

    # Update keys vector - remove both keys
    filter!(k -> k != p.ki && k != p.kj, p.keys_vec)

    # remove keys associated with old bond
    delete!(p.pos, p.ki)
    delete!(p.pos, p.kj)

    # If removing inside a pattern (type >= 2), extend pattern at tail
    if currentType >= 2 && p.patternTails !== nothing && haskey(p.patternTails, currentType)
        tailKey = p.patternTails[currentType]
        nextAfterTail = p.pos[tailKey][5]
        
        # Check if there are at least 2 base-type beads after the tail
        if nextAfterTail != 0 && p.pos[nextAfterTail][6] == Int32(1)
            nextNextAfterTail = p.pos[nextAfterTail][5]
            if nextNextAfterTail != 0 && p.pos[nextNextAfterTail][6] == Int32(1)
                # Swap two base-type beads to pattern type
                p.pos[nextAfterTail][6] = currentType
                p.pos[nextNextAfterTail][6] = currentType
                
                # Update patternTails: new tail is two beads forward
                p.patternTails[currentType] = nextNextAfterTail
            end
            # If not enough base-type beads, the move should have been rejected earlier
        end
    end

    return p
end

# Finds the bead coordinates, neighbording beads, and vector of the bond 
# to which the bead selected by p.moveIdx belongs
function identifyBond!(p::Polymer)

    @assert p.moveIdx >= 1 && p.moveIdx <= p.length begin
        "p.moveIdx must be between 1 and Species length (# of monomers)"
    end

    # keys_array = collect(keys(p.pos))  # TODO: this wastes allocs, fix it

    # pointers to the two beads (i and j) to be moved
    p.ki = p.keys_vec[p.moveIdx]
    p.kj = p.pos[p.ki][5]

    # previous and next beads along the chain
    p.kprev = p.pos[p.ki][4]
    p.knext = p.kj == 0 ? 0 : p.pos[p.kj][5]

    # bond edge coordinates
    p.bi = view(p.pos[p.ki], range(1, 3))
    p.bj = p.kj == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.kj], range(1, 3))

    # coordinates of previous and next beads relative to bond
    p.bprev = p.kprev == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.kprev], range(1, 3))
    p.bnext = p.knext == 0 ? view(SENTINEL, range(1, 3)) : view(p.pos[p.knext], range(1, 3))

    for i in range(1, 3)
        p.bondVec[i] = (p.bj[i] - p.bi[i]) % p.boxlength
    end

    return p
end

# Calls the appropriate operation (add, kink or remove) to execute a move
function makeMove!(p::Polymer)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end

    # bond addition
    if p.moveType == 1
        addBond!(p)

    # kink
    elseif p.moveType == 2
        kinkBond!(p)

    # bond removal  
    elseif p.moveType == 3 
        removeBond!(p)
    end

    return p
end

# Check if a pattern removal move is valid (has enough base-type beads to compensate)
function canRemoveFromPattern(p::Polymer, currentType::Int32)
    if currentType < 2 || p.patternTails === nothing || !haskey(p.patternTails, currentType)
        return true  # Not a pattern type, or no tracking - allow
    end
    
    tailKey = p.patternTails[currentType]
    nextAfterTail = p.pos[tailKey][5]
    
    if nextAfterTail == 0 || p.pos[nextAfterTail][6] != Int32(1)
        return false
    end
    
    nextNextAfterTail = p.pos[nextAfterTail][5]
    if nextNextAfterTail == 0 || p.pos[nextNextAfterTail][6] != Int32(1)
        return false
    end
    
    return true
end

# executes a move on polymermer s inside lattice l at position p.moveIdx, of type
# moveType (add, remove or kink) sampled according to the metropolis algorithm
function move!(p::Polymer)

    # sample a bond at random along the chain
    p.moveIdx = rand(1:p.length - 1)
    
    # find coordinates and keys of selected bond and adjacent beads
    identifyBond!(p)
    
    # don't move if last bead was selected
    if p.kj != 0
        
        # reject if bond spans type boundary (different monomer types)
        if p.pos[p.ki][6] != p.pos[p.kj][6]
            return false
        end
        
        # sample a move according to proposition probabilities
        proposeMove!(p)
        
        # For pattern removal moves, check if there are enough base-type beads
        currentType = p.pos[p.ki][6]
        if p.moveType == 3 && !canRemoveFromPattern(p, currentType)
            return false
        end
        
        # check for max length constraint and overlap
        if !overlap(p) && !(p.moveType == 1 && p.length == p.maxlength)
            
            # accept move with metropolis criterion
            if acceptMove(p)
                
                makeMove!(p)
                p.length = length(keys(p.pos))   # TODO: refactor for speed
                return true
            end
        end
    end

    return false
end

# delegates moving to parent polymer
move!(m::PolyMonomer) = move!(m.poly)


############# BulkMonomer move functions #############
# MARK: -BulkMonomer

function overlap(b::AbstractBulkMonomer)
    return b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] > 0
end

function identifyBond!(b::BulkMonomer)
    b.bi = view(b.pos[b.moveIdx], range(1, 3))

    return b
end

# Translates a particle to a new position set by b.newPos
function makeMove!(b::BulkMonomer)
    
    # update lattice sites
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = 0
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # update current position
    for i in range(1, 3)
        b.pos[b.moveIdx][i] = b.newPos[i]
    end

    @assert b.latS[b.pos[b.moveIdx][1], 
                   b.pos[b.moveIdx][2], 
                   b.pos[b.moveIdx][3]] == b.moveIdx begin
        "Invalid index-position correspondence"
    end

    return b
end

# Executes a move on a species of bulkMonomer monomers diluted in the box
function move!(b::BulkMonomer)

    # sample a particle at random
    b.moveIdx = rand(1:b.length)

    # find coordinates and keys of selected bond
    identifyBond!(b)
    
    # sample a move uniformly in each direction
    proposeMove!(b)

    if !overlap(b)
        if acceptMove(b)
            makeMove!(b)
            return true
        end
    end

    return false
end


############# FastBulkMonomer move functions #############
# MARK: -FastBulkMonomer

function identifyBond!(b::FastBulkMonomer)
    b.moveIdx = b.latS[b.bi[1], b.bi[2], b.bi[3]]
    b.moveIdxNew = b.latS[b.newPos[1], b.newPos[2], b.newPos[3]]

    return b
end

# Translates a particle to a new position set by b.newPos
function makeMove!(b::FastBulkMonomer)
    
    # update lattice sites
    b.latS[b.bi[1], b.bi[2], b.bi[3]] = b.moveIdxNew
    b.latS[b.newPos[1], b.newPos[2], b.newPos[3]] = b.moveIdx

    # swap old with new position
    if b.moveIdx != 0
        for i in range(1, 3)
            b.pos[b.moveIdx][i] = b.newPos[i]
        end

        @assert b.latS[b.pos[b.moveIdx][1], 
                       b.pos[b.moveIdx][2], 
                       b.pos[b.moveIdx][3]] == b.moveIdx begin
            "Invalid index-position correspondence"
        end
    end

    # swap new with old position
    if b.moveIdxNew != 0
        for i in range(1, 3)
            b.pos[b.moveIdxNew][i] = b.bi[i]
        end

        @assert b.latS[b.pos[b.moveIdxNew][1], 
                       b.pos[b.moveIdxNew][2], 
                       b.pos[b.moveIdxNew][3]] == b.moveIdxNew begin
            "Invalid index-position correspondence"
        end
    end

    return b
end

# Executes a move on a species of bulkMonomer monomers diluted in the box
function move!(b::FastBulkMonomer)

    # pick two sites on the lattice uniformly without replacement
    proposeMove!(b)
    
    if xor(b.moveIdx == 0, b.moveIdxNew == 0) && acceptMove(b)
        makeMove!(b)
        return true
    end

    return false
end