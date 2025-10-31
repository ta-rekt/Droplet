const PKINK = 0.8

# initializes bead locations of a long polymer in a doubly-linked list
export initPoly!
function initPoly!(s::Species)
    s.name = string("poly", s.id)
    s.pos = OrderedDict{Int32, Vector{Int32}}()
    v = zeros(Int32, 3)

    # line from center of one face to the other along x axis
    for i in 1:s.length
        v[1] = s.initPos === nothing ? i : i + s.initPos[1]
        v[2] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[2]
        v[3] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[3]
        s.pos[i] = [v..., i-1, i+1]  # last two indices are pointers to prev and next bead
    end

    # head and tail beads
    s.head = 1
    s.pos[s.head][4] = 0
    s.pos[s.length][5] = 0

    # key vector
    s.keys_vec = collect(keys(s.pos))

    s.tempVec = zeros(Int32, 3)
    s.bondVec = zeros(Int32, 3)
    s.moveVec = zeros(Int32, 3)
    s.newPosi = zeros(Int32, 3)
    s.newPosj = zeros(Int32, 3)
    s.moveTypes = zeros(Int32, 4)
    s.moveTypeProbs = [0,PKINK,0]
    s.nTypes = zeros(Int32, 3)

    return s
end

# initialize hamiltonian
export initHamPoly!
function initHamPoly!(s::Species)
    hamPoly!(s)

    return s
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initDilute!
function initDilute!(s::Species)

    s.name = string("dilute", s.id)
    s.pos = OrderedDict{Int32, Vector{Int32}}()
    s.head = 1

    for i in 1:s.length
        s.pos[i] = rand(1:s.boxlength, 5)
        s.pos[i][4] = i-1
        s.pos[i][5] = i+1
    end

    s.pos[1][4] = 0
    s.pos[s.length][5] = 0

    return s
end

# initializes lattice sites of each species from their bead locations
export initLattice!
function initLattice!(l::Lattice, ss::Vector{Species})

    l.nspecies = length(ss)   # update number of species
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int, Int}()

    # fill out lattice sites for each species
    for (i, s) in enumerate(ss)
        l.speciesDir[s.id] = i

        # give species s a view on its lattice sites
        s.latS = view(l.lat, i, :, :, :)

        for j in 1:s.length
            v = view(s.pos[j], 1:3)
            s.latS[v...] = 1
        end
    end

    return l
end