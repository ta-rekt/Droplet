############# Species initialization functions #############

######## Generic #######

# initialize hamiltonian
export initHam!
function initHam!(s::AbstractSpecies)
    ham!(s)
    return s
end


# initializes lattice sites from a list of presets
export initShape!
function initShape!(s::AbstractBulk)
    if s.initShape == "uniform"
        # sample uniformly randomly
        pos_idxs = sample(0:(s.boxlength^3-1), s.length, replace=false)
        
        for i in 1:s.length
            s.pos[i] = zeros(5)
            mx = mod(pos_idxs[i], s.boxlength^2)
            
            s.pos[i][1] = div(pos_idxs[i], s.boxlength^2) + 1 # x coordinate
            s.pos[i][2] = div(mx, s.boxlength) + 1            # y coordinate
            s.pos[i][3] = mod(mx, s.boxlength) + 1            # z coordinate
            s.pos[i][4] = i - 1
            s.pos[i][5] = i + 1
        end
        
    elseif s.initShape == "line"
        # line from initial position
        v = zeros(Int32, 3)
        
        for i in 1:s.length
            s.pos[i][1] = s.initPos === nothing ? i : i + s.initPos[1] - 1
            s.pos[i][2] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[2]
            s.pos[i][3] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[3]
            s.pos[i][4] = i - 1
            s.pos[i][5] = i + 1
        end
        
    elseif s.initShape == "box"
        # filled box shape
        idx = 1
        box_side = ceil(Int, s.length^(1/3))
        centerxyz = div(s.boxlength, 2)
        startPos = s.initPos === nothing ? [centerxyz, centerxyz, centerxyz] : s.initPos
        
        for x in 0:(box_side-1)
            for y in 0:(box_side-1)
                for z in 0:(box_side-1)
                    if idx > s.length
                        break
                    end
                    s.pos[idx] = zeros(5)
                    s.pos[idx][1] = startPos[1] + x - div(box_side, 2)
                    s.pos[idx][2] = startPos[2] + y - div(box_side, 2)
                    s.pos[idx][3] = startPos[3] + z - div(box_side, 2)
                    s.pos[idx][4] = idx - 1
                    s.pos[idx][5] = idx + 1
                    idx += 1
                end
                if idx > s.length
                    break
                end
            end
            if idx > s.length
                break
            end
        end
    end
    
    # head and tail beads
    s.pos[1][4] = 0
    s.pos[s.length][5] = 0
    
    return s
end

######## Polymer #######

# initializes bead locations of a long polymer in a doubly-linked list
export initSpecies!
function initSpecies!(p::Poly)
    p.name = string("poly", p.id)
    p.pos = OrderedDict{Int32, Vector{Int32}}()
    v = zeros(Int32, 3)

    # line from center of one face to the other along x axis
    for i in 1:p.length
        v[1] = p.initPos === nothing ? i : i + p.initPos[1]
        v[2] = p.initPos === nothing ? div(p.boxlength, 2) : p.initPos[2]
        v[3] = p.initPos === nothing ? div(p.boxlength, 2) : p.initPos[3]
        p.pos[i] = [v..., i-1, i+1]  # last two indices are pointers to prev and next bead
    end

    # head and tail beads
    p.head = 1
    p.pos[p.head][4] = 0
    p.pos[p.length][5] = 0

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

    return p
end


######## Bulk #######

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initSpecies!
function initSpecies!(b::Bulk)

    b.name = string("bulk", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()

    pos_idxs = sample(0:(b.boxlength^3-1), b.length, replace=false)

    for i in 1:b.length
        b.pos[i] = zeros(5)
        mx = mod(pos_idxs[i], b.boxlength^2)

        b.pos[i][1] = div(pos_idxs[i], b.boxlength^2) + 1 # x coordinate
        b.pos[i][2] = div(mx, b.boxlength) + 1            # y coordinate
        b.pos[i][3] = mod(mx, b.boxlength) + 1            # z coordinate
        b.pos[i][4] = i - 1
        b.pos[i][5] = i + 1
    end

    # head and tail beads
    b.pos[1][4] = 0
    b.pos[b.length][5] = 0

    b.nnsOld = 0
    b.nnsNew = 0
    b.newPos = zeros(Int32, 3)
    b.tempVec = zeros(Int32, 3)

    return b
end


######## FastBulk #######

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initSpecies!
function initSpecies!(b::FastBulk)

    b.name = string("bulk", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()

    pos_idxs = sample(0:(b.boxlength^3-1), b.length, replace=false)

    for i in 1:b.length
        b.pos[i] = zeros(5)
        mx = mod(pos_idxs[i], b.boxlength^2)

        b.pos[i][1] = div(pos_idxs[i], b.boxlength^2) + 1 # x coordinate
        b.pos[i][2] = div(mx, b.boxlength) + 1            # y coordinate
        b.pos[i][3] = mod(mx, b.boxlength) + 1            # z coordinate
        b.pos[i][4] = i - 1
        b.pos[i][5] = i + 1
    end

    # build set of movable particles
    b.movableIdxs = Set{Int32}()
    b.freeNNs = Set{Int32}()
    b.hastingsFactor = 1

    for i in 1:b.length
        if nn(b, b.pos[i][range(1,3)]) < 6
            push!(b.movableIdxs, i)
        end
    end

    # head and tail beads
    b.pos[1][4] = 0
    b.pos[b.length][5] = 0

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
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int, Int}()

    # fill out lattice sites for each species
    for (i, s) in enumerate(slist)
        
        s.boxlength = l.boxlength
        initSpecies!(s)

        l.speciesDir[s.id] = i

        # give species s a view on its lattice sites
        s.latS = view(l.lat, i, :, :, :)

        # TODO: apply this to all species
        if s <: AbstractFastBulk
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = j
            end
        else
            for j in 1:s.length
                s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
            end
        end

        # requires s.latS to be set
        if hasproperty(s, :ham)
            initHam!(s)
        end
    end

    return l
end