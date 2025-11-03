############# Polymer initialization functions #############

# initializes bead locations of a long polymer in a doubly-linked list
export init!
function init!(p::Poly)
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

# initialize hamiltonian
export initHam!
function initHam!(p::Poly)
    ham!(p)
    return p
end


############# Bulk initialization functions #############

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export init!
function init!(b::Bulk)

    b.name = string("bulk", b.id)
    b.pos = OrderedDict{Int32, Vector{Int32}}()
    
    for i in 1:b.length
        b.pos[i] = rand(1:b.boxlength, 5)
        b.pos[i][4] = i-1
        b.pos[i][5] = i+1
    end
    
    # head and tail beads
    b.pos[1][4] = 0
    b.pos[b.length][5] = 0

    return b
end


############# Lattice initialization functions #############

# initializes lattice sites of each species from their bead locations
export initLattice!
function initLattice!(l::Lattice, slist::Vector{<:AbstractSpecies})

    l.nspecies = length(slist)   # update number of species
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int, Int}()

    # fill out lattice sites for each species
    for (i, s) in enumerate(slist)
        l.speciesDir[s.id] = i

        # give species s a view on its lattice sites
        s.latS = view(l.lat, i, :, :, :)

        for j in 1:s.length
            s.latS[s.pos[j][1], s.pos[j][2], s.pos[j][3]] = 1
        end
    end

    return l
end