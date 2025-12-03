############# Species initialization functions #############

######## Generic #######

# initialize hamiltonian
export initHam!
function initHam!(s::AbstractSpecies)
    ham!(s)
    return s
end

# initial configurations
export initUniform!
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

export initLine!
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

export initBox!
function initBox!(pos::OrderedDict{Int32, Vector{Int32}}, 
                  latticelength::Int32,
                  length::Int32,
                  startPos::Union{Nothing, Vector{Int32}},
                  newIndex::Int = 1)

    # filled box shape
    idx = newIndex
    box_side = ceil(Int, length^(1/3))
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

# populate the set of movable particles of a species
export initMovableSet!
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

# initializes empty dict of bead locations and other containers
export initSpecies!
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

# initializes species configuration from a list of presets
export initShape!
function initShape!(s::AbstractPoly)
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

# initializes bead locations of a long polymer in a doubly-linked list
export initSpecies!
function initSpecies!(p::Poly)
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


######## Bulk #######

# initializes species configuration from a list of presets
export initShape!
function initShape!(s::AbstractBulk)
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
export initSpecies!
function initSpecies!(b::Bulk)

    b.name = string("bulk", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()
    
    initShape!(b)
    
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0
    b.nnsOld = 0
    b.nnsNew = 0
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.nnIdxs = zeros(Int32, 6)

    return b
end


######## FastBulk #######

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initSpecies!
function initSpecies!(b::FastBulk)

    b.name = string("bulk", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()

    initShape!(b)

    # build set of movable particles
    b.movableIdxs = Set{Int32}()
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1.0

    b.nnsOld = 0
    b.nnsNew = 0
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)
    b.nnIdxs = zeros(Int32, 6)

    return b
end


############# System initialization functions #############

# initializes lattice sites of each species from their bead locations
export init!
function init!(l::Lattice, slist::Vector{<:AbstractSpecies})

    l.nspecies = length(slist)   # update number of species
    l.lat = zeros((l.nspecies + 1, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int, Int}()

    # initialize ghost species
    gs = Species(;)
    gs.latS = view(l.lat, 1, :, :, :)
    gs.boxlength = l.boxlength

    initSpecies!(gs)

    # initialize each physical species
    for (i, s) in enumerate(slist)
        
        # give species s a view on its lattice sites, and on the ghost species
        s.latS = view(l.lat, i + 1, :, :, :)
        s.ghost = gs

        s.boxlength = l.boxlength
        initSpecies!(s)

        l.speciesDir[s.id] = i

        # TODO: apply this to all species
        if typeof(s) <: AbstractBulk
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = j
            end
        else
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
            end
        end

        if typeof(s) <: FastBulk
            initMovableSet!(s)
        end
    end

    for s in slist
        # requires s.latS to be set
        if hasproperty(s, :ham)
            initHam!(s)
        end
    end

    return l
end
