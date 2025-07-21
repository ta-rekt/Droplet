# initializes bead locations of a long polymer
export initPoly!
function initPoly!(s::Species)
    s.pos = zeros(Int32, s.length, 3)
    s.name = string("poly", s.id)

    # line from center of one face to the other along x axis
    for i in 1:s.length
        v = view(s.pos, i, :)
        v[1] = i
        v[2] = s.init_pos === nothing ? div(s.length, 2) : s.init_pos[2]
        v[3] = s.init_pos === nothing ? div(s.length, 2) : s.init_pos[3]
    end

    return s
end

# initializes beads positioned uniformly randomly (dilute and well-mixed)
export initDilute!
function initDilute!(s::Species)
    s.pos = rand(1:s.boxlength, s.length, 3)

    return s
end

# initializes lattice sites of each species from their bead locations
export initLattice!
function initLattice!(l::Lattice, species::Vector{Species})
    l.nspecies = max(l.nspecies, length(species))   # update number of species
    l.lat = zeros((l.nspecies, l.boxlength, l.boxlength, l.boxlength))

    # fill out lattice sites for each species
    for s in species
        latS = view(l.lat, s.id, :, :, :)   # assumes species IDs are unique, increment by 1 and start at 1 -> TEST

        for j in 1:s.length
            v = view(s.pos, j, :)
            latS[CartesianIndex((v[1], v[2], v[3]))] = 1
        end
    end

    return l
end