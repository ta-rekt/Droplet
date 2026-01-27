
export initHam!, initUniform!, initLine!, initBox!, initMovableSet!, 
       initMonomers!, initShape!, init!

############# Species initialization functions #############

######## Generic #######
# MARK: -Generic

# initialize hamiltonian
function initHam!(s::AbstractSpecies)
    ham!(s)
    return s
end

# initial configurations
function initUniform!(pos::OrderedDict{Int32, Vector{Int32}}, boxlength::Int32, length::Int32)
    # sample uniformly randomly
    pos_idxs = sample(0:(boxlength^3-1), length, replace=false)
    
    for i in 1:length
        pos[i] = zeros(5)
        mx = mod(pos_idxs[i], boxlength^2)
        
        pos[i][1] = div(pos_idxs[i], boxlength^2) + 1 # x coordinate
        pos[i][2] = div(mx, boxlength) + 1            # y coordinate
        pos[i][3] = mod(mx, boxlength) + 1            # z coordinate
        pos[i][4] = i - 1
        pos[i][5] = i + 1
    end

    return pos
end

function initLine!(pos::OrderedDict{Int32, Vector{Int32}}, 
                   initPos::Union{Nothing, Vector{Int32}}, 
                   boxlength::Int32,
                   length::Int32)

    # line from center of one face to the other along x axis
    for i in 1:length
        pos[i] = zeros(5)
        pos[i][1] = initPos === nothing ? i : i + initPos[1] - 1
        pos[i][2] = initPos === nothing ? div(boxlength, 2) : initPos[2]
        pos[i][3] = initPos === nothing ? div(boxlength, 2) : initPos[3]
        pos[i][4] = i - 1
        pos[i][5] = i + 1
    end

    return pos
end

function initBox!(pos::OrderedDict{Int32, Vector{Int32}}, 
                  latticelength::Int32,
                  length::Int32,
                  startPos::Union{Nothing, Vector{Int32}},
                  newIndex::Int32 = Int32(1))

    # filled box shape
    idx = newIndex
    box_side = ceil(Int32, length^(1/3))
    offset = Int32(div(box_side, 2))
    clippedPos = zeros(Int32, 3)
    rawPos = zeros(Int32, 3)
    
    for x in 0:(box_side-1)
        for y in 0:(box_side-1)
            for z in 0:(box_side-1)
                if idx >= length + newIndex
                    break
                end

                pos[idx] = zeros(5)
                rawPos = [startPos[1] + Int32(x) - offset,
                          startPos[2] + Int32(y) - offset,
                          startPos[3] + Int32(z) - offset]
                boxmod!(rawPos, clippedPos, Int32(1), latticelength)

                pos[idx][1] = clippedPos[1]
                pos[idx][2] = clippedPos[2]
                pos[idx][3] = clippedPos[3]
                pos[idx][4] = idx - 1
                pos[idx][5] = idx + 1
                idx += 1
            end
            if idx >= length + newIndex
                break
            end
        end
        if idx >= length + newIndex
            break
        end
    end
end

# Convenience method that converts Int64 to Int32
function initBox!(pos::OrderedDict{Int32, Vector{Int32}}, 
                  latticelength::Int32,
                  length::Int32,
                  startPos::Union{Nothing, Vector{Int32}},
                  newIndex::Int64)
    initBox!(pos, latticelength, length, startPos, Int32(newIndex))
end

# populate the set of movable particles of a species
function initMovableSet!(s::AbstractSpecies)
    for i in collect(keys(s.pos))
        nns = nn(s, s.pos[i][range(1,3)])
        if nns < 6
            push!(s.movableIdxs, i)
        end
    end

    return s
end

######## Species #######
# MARK: -Species

# initializes empty dict of bead locations and other containers
function initSpecies!(s::Species)
    s.name = string("species", s.id)
    s.pos = OrderedDict{Int32, Vector{Int32}}()

    # head and tail beads
    s.head = 1
    s.tempVec = zeros(Int32, 3)
    s.nnIdxs = zeros(Int32, 6)
    s.freeNNs = Set{Int32}()

    return s
end


######## Polymer #######
# MARK: -Polymer

"""
    initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)

Initialize subspecies for a polymer with bead-dependent interactions.

# Arguments
- `p::Polymer`: The polymer to initialize subspecies for
- `defaultChemPot::Float64`: Default chemical potential for all monomers (default: 0.0)
- `defaultIntPot::Float64`: Default interaction potential for all monomers (default: 0.0)

# Pattern options (set via p.pattern field):
- `:uniform`: All beads belong to subspecies 1 (default, backward compatible)
- `:alternating`: Alternates between subspecies (1, 2, 1, 2, ...)
- `:block`: Divides polymer into n equal blocks
- `:random`: Random assignment
- `:custom`: Use p.customPattern vector to assign beads

# Notes
- Each subspecies is initialized with the provided default potentials.
- Use intPot!(p, subspeciesIdx, value) to set subspecies-specific potentials.
- All subspecies share a synchronized intSpecies list and intMat (like AbstractBulk).
- When using :custom p.pattern, indices should be 1-based subspecies indices.
"""
function initMonomers!(p::Polymer; defaultChemPot::Float64=0.0, defaultIntPot::Float64=0.0)
    
    # Create subspecies with default values
    p.monomers = Vector{PolyMonomer}(undef, p.nMonomerTypes)
    for i in 1:p.nMonomerTypes
        p.monomers[i] = PolyMonomer(
            name = string(p.name, "_monomer", i),
            poly = p,
            chemPot = defaultChemPot,
            intPot = defaultIntPot
        )
    end
    
    # Initialize all subspecies with synced interaction structures
    # This ensures they all share the same intSpecies and intMat
    for i in 1:p.nMonomerTypes
        _ensureIntStructures!(p.monomers[i])
    end
    
    # Sync all subspecies to share the same matrix
    if p.nMonomerTypes > 1
        for i in 2:p.nMonomerTypes
            if p.monomers[1].intMat !== p.monomers[i].intMat
                _syncMatrices!(p.monomers[1], p.monomers[i])
            end
        end
    end
    
    # Initialize monomerIndex array for maxlength beads
    p.monomerIndex = zeros(Int32, p.maxlength)
    
    # Assign beads to subspecies based on p.pattern
    if p.pattern == :uniform
        # All beads are subspecies 1
        fill!(p.monomerIndex, 1)
    
    elseif p.pattern == :alternating
        # Alternate between subspecies
        for i in 1:p.maxlength
            p.monomerIndex[i] = mod1(i, p.nMonomerTypes)
        end
        
    elseif p.pattern == :block
        # Divide into n equal blocks
        blockSize = div(p.maxlength, p.nMonomerTypes)
        for i in 1:p.maxlength
            blockIdx = min(div(i - 1, blockSize) + 1, p.nMonomerTypes)
            p.monomerIndex[i] = blockIdx
        end
        
    elseif p.pattern == :random
        # Random assignment
        for i in 1:p.maxlength
            p.monomerIndex[i] = Int32(rand(1:p.nMonomerTypes))
        end
        
    elseif p.pattern == :custom
        if p.customPattern === nothing
            error("pattern=:custom requires p.customPattern vector to be provided")
        end
        if length(p.customPattern) != p.maxlength
            error("p.customPattern length ($(length(p.customPattern))) must equal p.maxlength ($(p.maxlength))")
        end
        # Validate all indices are in range
        for (i, idx) in enumerate(p.customPattern)
            if idx < 1 || idx > p.nMonomerTypes
                error("p.customPattern[$i] = $idx is out of range [1, $p.nMonomerTypes]")
            end
        end
        # Copy custom p.pattern
        p.monomerIndex .= p.customPattern
    else
        error("Unknown p.pattern: $pattern. Use :uniform, :alternating, :block, :random, or :custom")
    end
    
    return p
end

# initializes species configuration from a list of presets
function initShape!(s::Polymer)
    if s.initShape == "line" || s.initShape === nothing
        initLine!(s.pos, s.initPos, s.boxlength, s.length)
    else
        error("Invalid string for s.initShape")
    end
    
    # head and tail beads
    s.pos[1][4] = 0
    s.pos[s.length][5] = 0
    
    return s
end

# initializes bead locations of a long Polymer in a doubly-linked list
function initSpecies!(p::Polymer)
    p.name = string("poly", p.id)
    p.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(p)

    # head and tail beads
    p.head = 1

    # key vector
    p.keys_vec = collect(keys(p.pos))

    p.tempVec = zeros(Int32, 3)
    p.bondVec = zeros(Int32, 3)
    p.moveVec = zeros(Int32, 3)
    p.newPosi = zeros(Int32, 3)
    p.newPosj = zeros(Int32, 3)
    p.moveTypes = zeros(Int32, 4)
    p.moveTypeProbs = [0,PKINK,0]
    p.nTypes = zeros(Int32, 3)
    p.hastingsFactor = 1.0

    return p
end

######## BulkMonomer #######
# MARK: -BulkMonomer

# initializes species configuration from a list of presets
function initShape!(s::AbstractBulkMonomer)
    if s.initShape == "uniform" || s.initShape === nothing
        initUniform!(s.pos, s.boxlength, s.length)

    elseif s.initShape == "line"    
        initLine!(s.pos, s.initPos, s.boxlength, s.length)

    elseif s.initShape == "box"
        beadsPerBox = s.length
        totalBeads = beadsPerBox * s.initParams
        
        for i in 1:s.initParams
            initPos = zeros(Int32, 3)
            offset = div(s.boxlength, 2 * s.initParams)
            
            if !(s.initPos === nothing)
                initPos[1] = s.initPos[1] + offset * (2 * i - s.initParams - 1)
                initPos[2] = s.initPos[2]
                initPos[3] = s.initPos[3]
            else
                initPos[1] = div(s.boxlength, 2) + offset * (2 * i - s.initParams - 1)
                initPos[2] = div(s.boxlength, 2)
                initPos[3] = div(s.boxlength, 2)
            end

            startIdx = beadsPerBox * (i - 1) + 1
            initBox!(s.pos, s.boxlength, beadsPerBox, initPos, startIdx)
        end
        
        s.length = totalBeads   # total length = length of each box * nboxes
    else
        error("Invalid string for s.initShape")
    end
    
    # head and tail beads
    s.pos[1][4] = 0
    s.pos[s.length][5] = 0
    
    return s
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
function initSpecies!(b::BulkMonomer)

    b.name = string("BulkMonomer", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(b)

    _ensureIntStructures!(b)
    
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0
    b.nnsOld = 0
    b.nnsNew = 0
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.nnIdxs = zeros(Int32, 6)
    b.deltaH = 0.0
    b.deltaHP = 0.0

    if b.chemPot === nothing b.chemPot = 0.0 end
    if b.intPot === nothing b.intPot = 0.0 end

    return b
end


######## FastBulkMonomer #######
# MARK: -FastBulkMonomer

# initializes beads positioned uniformly randomly (dilute and well-mixed)
function initSpecies!(b::FastBulkMonomer)

    b.name = string("BulkMonomer", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()

    initShape!(b)

    _ensureIntStructures!(b)

    # build set of movable particles
    b.movableIdxs = Set{Int32}()
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0

    b.nnsOld = 0
    b.nnsNew = 0
    b.bi = zeros(Int32, 3)
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.moveVec = zeros(Int32, 2)
    b.nnIdxs = zeros(Int32, 6)
    b.deltaH = 0.0
    b.deltaHP = 0.0

    if b.chemPot === nothing b.chemPot = 0.0 end
    if b.intPot === nothing b.intPot = 0.0 end

    return b
end


############# System initialization functions #############
# MARK: -System

# initializes lattice sites of each species from their bead locations
function init!(l::Lattice, slist::Vector{<:AbstractSpecies})

    # Separate species that need lattice layers from monomers that delegate
    layeredSpecies = AbstractSpecies[]
    monomers = AbstractPolyMonomer[]
    polyParents = Set{AbstractPolymer}()
    
    for s in slist
        if s isa AbstractPolyMonomer
            push!(monomers, s)
            push!(polyParents, s.poly)
        else
            push!(layeredSpecies, s)
        end
    end
    
    # Combine regular species with parent polymers for lattice allocation
    allLayeredSpecies = vcat(layeredSpecies, collect(polyParents))
    
    @assert length(allLayeredSpecies) <= 26 "Too many species"

    l.nspecies = length(allLayeredSpecies)
    l.lat = zeros((l.nspecies + 1, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int32, Int32}()

    # initialize ghost species
    gs = Species(;)
    gs.latS = view(l.lat, 1, :, :, :)
    gs.boxlength = l.boxlength
    initSpecies!(gs)

    # initialize each layered species (regular species + parent polymers)
    for (i, s) in enumerate(allLayeredSpecies)
        
        # give species s a view on its lattice sites, and on the ghost species
        s.latS = view(l.lat, i + 1, :, :, :)
        s.ghost = gs
        s.boxlength = l.boxlength
        
        # Initialize species if not already done (parent polymers might be pre-initialized)
        if s.pos === nothing
            initSpecies!(s)
        end

        # give each species a unique letter suffix
        s.name = string(s.name, "_", ALPHABET[i])
        l.speciesDir[s.id] = i

        # Populate lattice
        if typeof(s) <: AbstractBulkMonomer
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = j
            end
        else
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
            end
        end

        if typeof(s) <: FastBulkMonomer
            initMovableSet!(s)
        end
    end

    # Initialize Hamiltonians for all layered species
    for s in allLayeredSpecies
        if hasproperty(s, :ham)
            initHam!(s)
        end
    end

    return l
end