# initializes bead locations of a long polymer in a doubly-linked list
export initPoly!
function initPoly!(s::Species)
    s.name = string("poly", s.id)
    v = zeros(Int32, 3)

    # line from center of one face to the other along x axis
    for i in 1:s.boxlength
        v[1] = i
        v[2] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[2]
        v[3] = s.initPos === nothing ? div(s.boxlength, 2) : s.initPos[3]
        s.pos[i] = [v..., i-1, i+1]  # last two indices are pointers to prev and next bead
    end

    s.pos[1][4] = 0
    s.pos[s.boxlength][5] = 0

    return s
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initDilute!
function initDilute!(s::Species)

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
function initLattice!(l::Lattice, species::Vector{Species})
    l.nspecies = max(l.nspecies, length(species))   # update number of species
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))
    l.speciesDir = Dict{Int, Int}()

    # fill out lattice sites for each species
    for (i, s) in enumerate(species)
        l.speciesDir[s.id] = i
        latS = view(l.lat, i, :, :, :) 

        for j in 1:s.length
            v = view(s.pos[j], 1:3)
            latS[v...] = 1
        end
    end

    return l
end